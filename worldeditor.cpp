#include "precomp.h"
#include "worldeditor.h"

#define MAX_ALLOCATION 512000000 // in bytes

WorldEditor::WorldEditor()
{
	tempBricks = (PAYLOAD*)_aligned_malloc(CHUNKCOUNT * CHUNKSIZE, 64);
	tempGrid = (uint*)_aligned_malloc(GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * 4, 64);
	tempBrickInfo = (BrickInfo*)_aligned_malloc(BRICKCOUNT * sizeof(BrickInfo), 64);
	stateHead = (State*)calloc(1, sizeof(State));
	stateCurrent = stateTail = stateHead;

	memset(tempBricks, 0, CHUNKCOUNT * CHUNKSIZE);
	memset(tempGrid, 0, GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));
	memset(tempBrickInfo, BRICKSIZE, BRICKCOUNT * sizeof(BrickInfo));

	loadedTiles.push_back(LoadTile("assets/flower.vox"));
	tileIdx = loadedTiles.front();
}

#pragma region InputHandling
void WorldEditor::MouseMove(int x, int y)
{
	if (mousePos.x == x && mousePos.y == y) return;

	// Keep track of mouse position so if re-enabled, the selected brick is accurate
	mousePos.x = x, mousePos.y = y;

	if (!enabled) return;

	// Just hovering
	if (gesture.state == GestureState::GESTURE_POSSIBLE)
	{
		// Reset the selected voxel aabb 
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();
		return;
	}

	if (gesture.mode & GestureMode::GESTURE_MULTI)
	{
		MultiAddRemove();
	}
	else 
	{
		auto oldBox = selectedBricks.box;
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();

		// Ignore when mouse hovers over the same brick to avoid multiple add/removes
		if (oldBox == selectedBricks.box) return;

		if (gesture.mode & GestureMode::GESTURE_REMOVE) RemoveSelectedBrick();
		else AddSelectedBrick();
	}
	UpdateSelectedBrick();
}

void WorldEditor::KeyDown(int key)
{
	if (!enabled) return;

	switch (key) {
	case GLFW_KEY_LEFT_CONTROL:
		selectedKeys |= GestureKey::GESTURE_CTRL;
		UpdateGestureMode();
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();
		break;
	case GLFW_KEY_LEFT_SHIFT:
		selectedKeys |= GestureKey::GESTURE_SHIFT;
		UpdateGestureMode();
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();
		break;
	case GLFW_KEY_Z:
		if (selectedKeys == GestureKey::GESTURE_CTRL)
		{
			Undo();
			selectedBricks.box = aabb{};
			UpdateSelectedBrick();
		}
		break;
	case GLFW_KEY_Y:
		if (selectedKeys == GestureKey::GESTURE_CTRL)
		{
			Redo();
			selectedBricks.box = aabb{};
			UpdateSelectedBrick();
		}
		break;
	default:
		break;
	}
}

void WorldEditor::KeyUp(int key)
{
	if (!enabled) return;

	switch (key) {
	case GLFW_KEY_LEFT_CONTROL:
		selectedKeys ^= GestureKey::GESTURE_CTRL;
		UpdateGestureMode();
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();
		break;
	case GLFW_KEY_LEFT_SHIFT:
		selectedKeys ^= GestureKey::GESTURE_SHIFT;
		UpdateGestureMode();
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();
		break;
	default:
		break;
	}
}

void WorldEditor::MouseDown(int mouseButton)
{
	if (!enabled) return;

	if (gesture.state != GestureState::GESTURE_POSSIBLE) return;

	switch (mouseButton) {
	case GLFW_MOUSE_BUTTON_LEFT:
		selectedButtons |= GestureButton::GESTURE_LMB;
		break;
	case GLFW_MOUSE_BUTTON_RIGHT:
		selectedButtons |= GestureButton::GESTURE_RMB;
		break;
	default:
		break;
	}

	UpdateGestureMode();

	if (gesture.buttons & GestureButton::GESTURE_LMB)
	{
		// Clear the list of updated bricks for the new gesture
		updatedBricks.clear();
		gesture.state = GestureState::GESTURE_ACTIVE;

		// Save a backup version of the grid/brick/info 
		World& world = *GetWorld();
		memcpy(tempBricks, world.GetBrick(), CHUNKCOUNT * CHUNKSIZE);
		memcpy(tempGrid, world.GetGrid(), GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));
		memcpy(tempBrickInfo, world.GetBrickInfo(), GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));

		// Add or remove the intial brick where the mouse is hovered
		if (gesture.mode & GestureMode::GESTURE_REMOVE) RemoveSelectedBrick();
		else AddSelectedBrick();

		// Set our anchor for multi add/remove
		selectedBricks.anchor = selectedBricks.box;
	}

	UpdateSelectedBrick();
}

void WorldEditor::MouseUp(int mouseButton)
{
	if (!enabled) return;

	bool gestureFinished = false;
	switch (mouseButton) {
	case GLFW_MOUSE_BUTTON_LEFT:
		if (gesture.buttons & GestureButton::GESTURE_LMB) gestureFinished = true;
		selectedButtons ^= GestureButton::GESTURE_LMB;
		break;
	case GLFW_MOUSE_BUTTON_RIGHT:
		selectedButtons ^= GestureButton::GESTURE_RMB;
		break;
	default:
		break;
	}

	// If we've completed a gesture, save the state for undo/redo purposes
	if (gestureFinished)
	{
		SaveState();
		gesture.state = GestureState::GESTURE_POSSIBLE;
		UpdateGestureMode();

		// Reset the selected voxel aabb 
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();
	}
}
#pragma endregion

#pragma region Editing
void WorldEditor::MultiAddRemove()
{
	aabb oldBox = selectedBricks.box;
	selectedBricks.box = selectedBricks.anchor;
	UpdateSelectedBrick();

	// Check to see if the selected bounding box has changed
	if (!(oldBox == selectedBricks.box))
	{
		World& world = *GetWorld();
		unsigned short* brick = world.GetBrick();
		uint* grid = world.GetGrid();

		// Calculate the symmetrical difference
		// Subtract the intersextion from oldbox is voxels to remove
		// Subtracting intersection from newbox is voxels to add
		aabb intersection = oldBox.Intersection(selectedBricks.box);

		// We've updated the selected bricks box so now it doesn't cover the same bricks
		// For all updated bricks, put back the old bricks values
		for (int bx = oldBox.bmax3.x; bx >= oldBox.bmin3.x; bx--)
			for (int by = oldBox.bmax3.y; by >= oldBox.bmin3.y; by--)
				for (int bz = oldBox.bmax3.z; bz >= oldBox.bmin3.z; bz--)
				{
					__m128 b4 = _mm_setr_ps(bx, by, bz, 0);
					if (intersection.Contains(b4))
						continue;

					const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
					const uint tempGridVal = tempGrid[cellIdx];
					const uint curGridVal = grid[cellIdx];

					// If we're reintroducing what was an empty brick, just remove the current one
					if ((tempGridVal & 1) == 0)
					{
						world.RemoveBrick(bx, by, bz);
						continue;
					}

					// Get current and new brick offsets (g1)
					uint tempBrickOffset = (tempGridVal >> 1);
					uint curBrickOffset = (curGridVal >> 1);

					uint brickIdx = tempBrickOffset;
					uint gridValue = tempGridVal;

					// If the current brick is empty we need to create a NewBrick
					if ((curGridVal & 1) == 0) brickIdx = world.NewBrick(), gridValue = (brickIdx << 1) | 1;

					// Copy the brick from the saved temporary brick buffer to our current brick buffer 
					memcpy(brick + brickIdx * BRICKSIZE, tempBricks + tempBrickOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
					grid[cellIdx] = gridValue;
					world.GetBrickInfo()[brickIdx].zeroes = tempBrickInfo[tempBrickOffset].zeroes;
					world.Mark(brickIdx);
				}

		updatedBricks.clear();
		// Either add or remove the bricks in the new selected bricks aabb
		for (int bx = selectedBricks.box.bmax3.x; bx >= selectedBricks.box.bmin3.x; bx--)
			for (int by = selectedBricks.box.bmax3.y; by >= selectedBricks.box.bmin3.y; by--)
				for (int bz = selectedBricks.box.bmax3.z; bz >= selectedBricks.box.bmin3.z; bz--)
				{
					updatedBricks.push_back((make_int3)(bx, by, bz));
					__m128 b4 = _mm_setr_ps(bx, by, bz, 0);
					if (intersection.Contains(b4))
						continue;

					if (gesture.mode & GestureMode::GESTURE_REMOVE)
					{
						world.RemoveBrick(bx, by, bz);
					}
					else
					{
						world.DrawTile(tileIdx, bx, by, bz);
					}
				}

	}
}

void WorldEditor::AddSelectedBrick()
{
	World& world = *GetWorld();
	for (int bx = selectedBricks.box.bmin[0]; bx <= selectedBricks.box.bmax[0]; bx++)
		for (int by = selectedBricks.box.bmin[1]; by <= selectedBricks.box.bmax[1]; by++)
			for (int bz = selectedBricks.box.bmin[2]; bz <= selectedBricks.box.bmax[2]; bz++)
			{
				// Avoid adding the same brick multiple times
				int3 brickPos = make_int3(bx, by, bz);
				if (std::find(updatedBricks.begin(), updatedBricks.end(), brickPos) == updatedBricks.end())
				{
					// Just draw the tile regardless for now
					// Potentially should make this return a boolean value to determine if the tile was already there
					updatedBricks.push_back(brickPos);
					world.DrawTile(tileIdx, bx, by, bz);
				}
			}
}

void WorldEditor::RemoveSelectedBrick()
{
	World& world = *GetWorld();
	for (int bx = selectedBricks.box.bmin3.x; bx <= selectedBricks.box.bmax3.x; bx++)
		for (int by = selectedBricks.box.bmin3.y; by <= selectedBricks.box.bmax3.y; by++)
			for (int bz = selectedBricks.box.bmin3.z; bz <= selectedBricks.box.bmax3.z; bz++)
			{
				// Avoid removing the same brick multiple times
				int3 brickPos = make_int3(bx, by, bz);
				if (std::find(updatedBricks.begin(), updatedBricks.end(), brickPos) == updatedBricks.end())
				{
					const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
					const uint curGridVal = world.GetGrid()[cellIdx];
					if ((curGridVal & 1) == 0) continue;

					updatedBricks.push_back(brickPos);
					world.RemoveBrick(bx, by, bz);
				}
			}
}
#pragma endregion

#pragma region State

// Reset the world editor so when re-enabled we're clean 
void WorldEditor::ResetEditor()
{
	gesture = Gesture{};
	selectedBricks = Selected{};
	selectedKeys = GestureKey::GESTURE_NO_KEYS;
	selectedButtons = GestureButton::GESTURE_NO_BUTTONS;
}

// For Undo we want to copy the state's old values into the current world
void WorldEditor::Undo()
{
	// At head of state linked list, no more undos
	if (!stateCurrent->prevState || gesture.state != GestureState::GESTURE_POSSIBLE || !undoEnabled) return;

	World& world = *GetWorld();
	uint* grid = world.GetGrid();
	PAYLOAD* bricks = world.GetBrick();
	BrickInfo* brickInfo = world.GetBrickInfo();

	int numBricks = stateCurrent->numBricks;
	for (int idx = 0; idx < numBricks; idx++)
	{
		int3 b = stateCurrent->updatedBricks[idx];
		const uint cellIdx = b.x + b.z * GRIDWIDTH + b.y * GRIDWIDTH * GRIDDEPTH;

		const uint oldGridVal = stateCurrent->oldGridVals[idx];
		const uint newGridVal = stateCurrent->newGridVals[idx];

		// If we're restoring what was an empty brick, just remove the current one
		if ((oldGridVal & 1) == 0)
		{
			world.RemoveBrick(b.x, b.y, b.z);
			continue;
		}

		// Get old and new brick offsets (g1)
		uint oldBrickOffset = (oldGridVal >> 1);
		uint newBrickOffset = (newGridVal >> 1);

		uint brickIdx = oldBrickOffset;
		uint gridValue = oldGridVal;

		// If the current brick is empty we need to create a NewBrick
		if ((newGridVal & 1) == 0) brickIdx = world.NewBrick(), gridValue = (brickIdx << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(bricks + (brickIdx * BRICKSIZE), stateCurrent->oldBricks + (idx * BRICKSIZE), BRICKSIZE * PAYLOADSIZE);
		grid[cellIdx] = gridValue;
		brickInfo[brickIdx].zeroes = stateCurrent->oldBrickZeroes[idx];
		world.Mark(brickIdx);
	}

	stateCurrent = stateCurrent->prevState;
}

// For Redo we want to copy the state's new values into the current world
void WorldEditor::Redo()
{
	// At tail of state linked list, no more undos
	if (!stateCurrent->nextState || gesture.state != GestureState::GESTURE_POSSIBLE || !undoEnabled) return;
	stateCurrent = stateCurrent->nextState;

	World& world = *GetWorld();
	uint* grid = world.GetGrid();
	PAYLOAD* bricks = world.GetBrick();
	BrickInfo* brickInfo = world.GetBrickInfo();

	int numBricks = stateCurrent->numBricks;
	for (int idx = 0; idx < numBricks; idx++)
	{
		int3 b = stateCurrent->updatedBricks[idx];
		const uint cellIdx = b.x + b.z * GRIDWIDTH + b.y * GRIDWIDTH * GRIDDEPTH;

		const uint oldGridVal = stateCurrent->oldGridVals[idx];
		const uint newGridVal = stateCurrent->newGridVals[idx];

		// If we're restoring what was an empty brick, just remove the current one
		if ((newGridVal & 1) == 0)
		{
			world.RemoveBrick(b.x, b.y, b.z);
			continue;
		}

		// Get old and new brick offsets (g1)
		uint oldBrickOffset = (oldGridVal >> 1);
		uint newBrickOffset = (newGridVal >> 1);

		uint brickIdx = newBrickOffset;
		uint gridValue = newGridVal;

		// If the current brick is empty we need to create a NewBrick
		if ((oldGridVal & 1) == 0) brickIdx = world.NewBrick(), gridValue = (brickIdx << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(bricks + (brickIdx * BRICKSIZE), stateCurrent->newBricks + (idx * BRICKSIZE), BRICKSIZE * PAYLOADSIZE);
		grid[cellIdx] = gridValue;
		brickInfo[brickIdx].zeroes = stateCurrent->oldBrickZeroes[idx];
		world.Mark(brickIdx);
	}

}

// Called after a gesture has been completed so we can update the history
void WorldEditor::SaveState()
{
	World& world = *GetWorld();
	uint* grid = world.GetGrid();
	PAYLOAD* bricks = world.GetBrick();
	BrickInfo* brickInfo = world.GetBrickInfo();

	if (!CreateNewState())
	{
		undoEnabled = false;
		return;
	}

	int numBricks = updatedBricks.size();
	stateCurrent->numBricks = numBricks;
	stateCurrent->newBricks = (PAYLOAD*)_aligned_malloc(numBricks * BRICKSIZE * PAYLOADSIZE, 64);
	stateCurrent->oldBricks = (PAYLOAD*)_aligned_malloc(numBricks * BRICKSIZE * PAYLOADSIZE, 64);
	stateCurrent->newGridVals = (uint*)_aligned_malloc(numBricks * 4, 64);
	stateCurrent->oldGridVals = (uint*)_aligned_malloc(numBricks * 4, 64);
	stateCurrent->newBrickZeroes = (uint*)_aligned_malloc(numBricks * sizeof(uint), 64);
	stateCurrent->oldBrickZeroes = (uint*)_aligned_malloc(numBricks * sizeof(uint), 64);
	stateCurrent->updatedBricks = (int3*)_aligned_malloc(numBricks * sizeof(int3), 64);

	for (int idx = 0; idx < numBricks; idx++)
	{
		int3 brick = updatedBricks[idx];
		const uint cellIdx = brick.x + brick.z * GRIDWIDTH + brick.y * GRIDWIDTH * GRIDDEPTH;

		const uint tempGridValue = tempGrid[cellIdx];
		const uint gridValue = grid[cellIdx];

		uint tempBrickOffset = (tempGridValue >> 1);
		uint brickOffset = (gridValue >> 1);

		// Old brick is empty
		if ((tempGridValue & 1) == 0)
		{
			memset(stateCurrent->oldBricks + (idx * BRICKSIZE), 0, BRICKSIZE * PAYLOADSIZE);
		}
		else
		{
			memcpy(stateCurrent->oldBricks + (idx * BRICKSIZE), tempBricks + tempBrickOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
		}
		
		// New brick is empty
		if ((gridValue & 1) == 0)
		{
			memset(stateCurrent->newBricks + (idx * BRICKSIZE), 0, BRICKSIZE * PAYLOADSIZE);
		}
		else
		{
			memcpy(stateCurrent->newBricks + (idx * BRICKSIZE), bricks + brickOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
		}

		stateCurrent->oldGridVals[idx] = tempGridValue;
		stateCurrent->newGridVals[idx] = gridValue;

		stateCurrent->oldBrickZeroes[idx] = tempBrickInfo[tempBrickOffset].zeroes;
		stateCurrent->newBrickZeroes[idx] = brickInfo[brickOffset].zeroes;

		stateCurrent->updatedBricks[idx] = brick;
	}


	// Add the amount of allocated memory 
	allocatedUndo += 2 * (numBricks * BRICKSIZE * PAYLOADSIZE);
	allocatedUndo += 2 * (numBricks * 4);
	allocatedUndo += 2 * (numBricks * sizeof(uint));
	allocatedUndo += (numBricks * sizeof(int3));


	// Check to see if the state has actually changed....if not don't add it to the history
	if (memcmp(stateCurrent->oldBricks, stateCurrent->newBricks, numBricks * BRICKSIZE * PAYLOADSIZE) == 0)
	{
		stateCurrent = stateCurrent->prevState;
		DeleteState(stateCurrent->nextState);
		stateCurrent->nextState = NULL;
		stateTail = stateCurrent;
		return;
	}

	CheckMemoryAllowance();
}

// Called at start of a gesture to establish a new state for the upcoming changes
bool WorldEditor::CreateNewState()
{
	State* newState = (State*)calloc(1, sizeof(State));
	if (!newState) return false;

	newState->prevState = stateCurrent;
	newState->nextState = NULL;

	// If we perform an action after many undos, free all undone states
	while (stateCurrent != stateTail)
	{
		stateTail = stateTail->prevState;
		DeleteState(stateTail->nextState);
		stateTail->nextState = NULL;
	}

	stateCurrent->nextState = newState;	
	stateCurrent = stateTail = newState;

	return true;
}

// Free up all allocated memory for our state
void WorldEditor::DeleteState(State* state)
{
	int numBricks = state->numBricks;

	_aligned_free(state->oldBricks);
	_aligned_free(state->newBricks);
	_aligned_free(state->oldGridVals);
	_aligned_free(state->newGridVals);
	_aligned_free(state->oldBrickZeroes);
	_aligned_free(state->newBrickZeroes);
	_aligned_free(state->updatedBricks);
	free(state);

	allocatedUndo -= 2 * (numBricks * BRICKSIZE * PAYLOADSIZE);
	allocatedUndo -= 2 * (numBricks * 4);
	allocatedUndo -= 2 * (numBricks * sizeof(uint));
	allocatedUndo -= (numBricks * sizeof(int3));
}
#pragma endregion

// Update which brick is currently selected by the mouse cursor 
void WorldEditor::UpdateSelectedBrick()
{
	World& world = *GetWorld();
	RenderParams& params = world.GetRenderParams();

	// setup primary ray for pixel [x,y] 
	float u = (float)mousePos.x / SCRWIDTH;
	float v = (float)mousePos.y / SCRHEIGHT;

	const float2 uv = make_float2(mousePos.x * params.oneOverRes.x, mousePos.y * params.oneOverRes.y);
	const float3 P = params.p0 + (params.p1 - params.p0) * uv.x + (params.p2 - params.p0) * uv.y;

	Ray ray;
	ray.O = params.E;
	ray.D = normalize(P - params.E);
	ray.t = 1e34f;

	// Trace the ray using the previous grid/brick state if gesture is active
	Intersection intersection = (gesture.state == GestureState::GESTURE_ACTIVE) ? Trace(ray, tempBricks, tempGrid) : Trace(ray);

	// Check to see if we hit the world grid
	if (intersection.GetVoxel() == 0)
	{
		float3 radius = make_float3(512, 512, 512);

		// Move to the box's reference frame. This is unavoidable and un-optimizable.
		float3 rayOrigin = ray.O - radius;


		// Winding direction: -1 if the ray starts inside of the box (i.e., and is leaving), +1 if it is starting outside of the box
		float3 res = fabs(rayOrigin) * (1.0 / radius);
		float max = res.x;
		if (res.y > max) max = res.y;
		if (res.z > max) max = res.z;
		float winding = (max < 1.0) ? -1.0 : 1.0;


		// We'll use the negated sign of the ray direction in several places, so precompute it.
		// The sign() instruction is fast...but surprisingly not so fast that storing the result
		// temporarily isn't an advantage.
		float3 sgn;
		sgn.x = signbit(ray.D.x) ? 1 : -1;
		sgn.y = signbit(ray.D.y) ? 1 : -1;
		sgn.z = signbit(ray.D.z) ? 1 : -1;

		// Ray-plane intersection. For each pair of planes, choose the one that is front-facing
		// to the ray and compute the distance to it.
		float3 distanceToPlane = radius * -1.0 * sgn - rayOrigin;
		distanceToPlane *= (1.0 / ray.D);

		// Perform all three ray-box tests and cast to 0 or 1 on each axis. 
		// Use a macro to eliminate the redundant code (no efficiency boost from doing so, of course!)
		// Could be written with 
		int3 test;


		test.x = (distanceToPlane.x >= 0.0) && LessThan(fabs(make_float2(rayOrigin.y, rayOrigin.z) + make_float2(ray.D.y, ray.D.z) * distanceToPlane.x), make_float2(radius.y, radius.z));
		test.y = (distanceToPlane.y >= 0.0) && LessThan(fabs(make_float2(rayOrigin.z, rayOrigin.x) + make_float2(ray.D.z, ray.D.x) * distanceToPlane.y), make_float2(radius.z, radius.x));
		test.z = (distanceToPlane.z >= 0.0) && LessThan(fabs(make_float2(rayOrigin.x, rayOrigin.y) + make_float2(ray.D.x, ray.D.y) * distanceToPlane.z), make_float2(radius.x, radius.y));

		// CMOV chain that guarantees exactly one element of sgn is preserved and that the value has the right sign
		sgn = test.x ? make_float3(sgn.x, 0.0, 0.0) : (test.y ? make_float3(0.0, sgn.y, 0.0) : make_float3(0.0, 0.0, test.z ? sgn.z : 0.0));

		// At most one element of sgn is non-zero now. That element carries the negative sign of the 
		// ray direction as well. Notice that we were able to drop storage of the test vector from registers,
		// because it will never be used again.

		// Mask the distance by the non-zero axis
		// Dot product is faster than this CMOV chain, but doesn't work when distanceToPlane contains nans or infs. 
		//
		float distance = (sgn.x != 0.0) ? distanceToPlane.x : ((sgn.y != 0.0) ? distanceToPlane.y : distanceToPlane.z);

		// Normal must face back along the ray. If you need
		// to know whether we're entering or leaving the box, 
		// then just look at the value of winding. If you need
		// texture coordinates, then use box.invDirection * hitPoint.
		float3 normal = sgn;

		bool hit = (sgn.x != 0) || (sgn.y != 0) || (sgn.z != 0);
		if (!hit) return;


		float3 hitpoint = ray.O + ray.D * distance;

		// Get position inside of the voxel to determine brick location
		float3 voxelPos = hitpoint + 0.5 * normal;

		float3 brickPos = make_float3((int)voxelPos.x / BRICKDIM, (int)voxelPos.y / BRICKDIM, (int)voxelPos.z / BRICKDIM);
		selectedBricks.box.Grow(brickPos);
	}
	else
	{
		float t = intersection.GetDistance();
		float3 N = intersection.GetNormal();
		float3 hitPoint = ray.O + ray.D * t;

		// Get position inside of the voxel to determine brick location
		float3 voxelPos = hitPoint - 0.1 * N;
		float3 brickPos = make_float3((int)voxelPos.x / BRICKDIM, (int)voxelPos.y / BRICKDIM, (int)voxelPos.z / BRICKDIM);

		if (gesture.mode & GestureMode::GESTURE_REMOVE)
		{
			selectedBricks.box.Grow(brickPos);
		}
		else
		{
			selectedBricks.box.Grow(brickPos + N);
		}
	}

	// Update rendering params to trace the selected bricks outline
	params.selectedMin = selectedBricks.box.bmin3 * BRICKDIM;
	params.selectedMax = selectedBricks.box.bmax3 * BRICKDIM + make_float3(BRICKDIM, BRICKDIM, BRICKDIM);
}

// Sets the current mode based on mouse button and selected keys
void WorldEditor::UpdateGestureMode()
{
	if (gesture.state != GestureState::GESTURE_POSSIBLE) return;
	gesture.buttons = selectedButtons;
	gesture.keys = selectedKeys;

	// Set the default state to add
	gesture.mode = GestureMode::GESTURE_ADD;
	if (gesture.keys & GestureKey::GESTURE_SHIFT) gesture.mode |= GestureMode::GESTURE_REMOVE;
	if (gesture.keys & GestureKey::GESTURE_CTRL) gesture.mode |= GestureMode::GESTURE_MULTI;
}

// Check to see if we've gone above the maximum limit of allocated memory for undo states
void WorldEditor::CheckMemoryAllowance()
{
	while (allocatedUndo > MAX_ALLOCATION && stateHead->nextState)
	{
		// Delete the state we just allocated
		if (stateHead->nextState == stateCurrent)
			stateCurrent = stateTail = stateHead;

		State* newState = stateHead->nextState->nextState;
		DeleteState(stateHead->nextState);
		stateHead->nextState = newState;
		if (newState) newState->prevState = stateHead;
	}
}