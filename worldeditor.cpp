#include "precomp.h"
#include "worldeditor.h"

WorldEditor::WorldEditor()
{
	tempBricks = (PAYLOAD*)_aligned_malloc(CHUNKCOUNT * CHUNKSIZE, 64);
	tempGrid = (uint*)_aligned_malloc(GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * 4, 64);
	memset(tempBricks, 0, CHUNKCOUNT * CHUNKSIZE);
	memset(tempGrid, 0, GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));
}

void WorldEditor::MouseMove(int x, int y)
{
	mousePos.x = x, mousePos.y = y;

	// just hovering
	if (gesture.state == GestureState::GESTURE_POSSIBLE)
	{
		// Reset the selected voxel aabb 
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();
		return;
	}

	if (gesture.state == GestureState::GESTURE_START)
	{
		gesture.state = GestureState::GESTURE_UPDATE;
	}

	// Default brush behaviour
	if (gesture.key == GestureKey::GESTURE_DEFAULT)
	{
		selectedBricks.box = aabb{};
		UpdateSelectedBrick();

		if (gesture.button == GestureButton::GESTURE_LMB)
		{
			if (gesture.mode == GestureMode::GESTURE_ADD)
			{
				AddBrick();
			}

			if (gesture.mode == GestureMode::GESTURE_SUBTRACT)
			{
				RemoveBrick();
			}
		}
	}
	else if (gesture.key == GestureKey::GESTURE_CTRL)
	{
		aabb oldBox = selectedBricks.box;
		selectedBricks.box = selectedBricks.anchor;
		UpdateSelectedBrick();

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
						const uint value = tempGrid[cellIdx];
						uint brickOffset = (value >> 1);

						memcpy(brick + brickOffset * BRICKSIZE, tempBricks + brickOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
						grid[cellIdx] = value;
						world.Mark(brickOffset);
					}

			// Either add or remove the bricks in the new selected bricks aabb
			for (int bx = selectedBricks.box.bmax3.x; bx >= selectedBricks.box.bmin3.x; bx--)
				for (int by = selectedBricks.box.bmax3.y; by >= selectedBricks.box.bmin3.y; by--)
					for (int bz = selectedBricks.box.bmax3.z; bz >= selectedBricks.box.bmin3.z; bz--)
					{
						__m128 b4 = _mm_setr_ps(bx, by, bz, 0);
						if (intersection.Contains(b4))
							continue;

						if (gesture.mode == GestureMode::GESTURE_ADD)
						{
							world.SetBrick(bx * BRICKDIM, by * BRICKDIM, bz * BRICKDIM, WHITE);

						}

						if (gesture.mode == GestureMode::GESTURE_SUBTRACT)
						{
							world.RemoveBrick(bx, by, bz);
						}
					}

		}
	}

}

void WorldEditor::KeyDown(int key)
{
	// Make sure ctrl button isn't pressed during another gesture
	if (gesture.state != GestureState::GESTURE_POSSIBLE)
	{
		return;
	}

	if (key == GLFW_KEY_LEFT_CONTROL)
	{
		gesture.key = GestureKey::GESTURE_CTRL;
	}

	if (key == GLFW_KEY_LEFT_SHIFT)
	{
		if (gesture.mode == GestureMode::GESTURE_ADD)
		{
			gesture.mode = GestureMode::GESTURE_SUBTRACT;
		}
		else if (gesture.mode == GestureMode::GESTURE_SUBTRACT)
		{
			gesture.mode = GestureMode::GESTURE_ADD;
		}
	}
}

void WorldEditor::KeyUp(int key)
{
	if (key == GLFW_KEY_LEFT_CONTROL)
	{
		// If the ctrl button wasn't pressed at beginning of gesture ignore the keyup
		if (!(gesture.key == GestureKey::GESTURE_CTRL))
		{
			return;
		}

		gesture.key = GestureKey::GESTURE_DEFAULT;
	}
}

void WorldEditor::MouseDown(int mouseButton)
{
	if (gesture.state != GestureState::GESTURE_POSSIBLE)
	{
		return;
	}

	gesture.state = GestureState::GESTURE_START;
	gesture.button = (GestureButton)mouseButton;

	// For now, add a brick to the face of the selected voxel the cursor is on
	if (gesture.button == GestureButton::GESTURE_LMB)
	{
		World& world = *GetWorld();
		memcpy(tempBricks, world.GetBrick(), CHUNKCOUNT * CHUNKSIZE);
		memcpy(tempGrid, world.GetGrid(), GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));

		if (gesture.mode == GestureMode::GESTURE_ADD)
		{
			AddBrick();

		}
		else if (gesture.mode == GestureMode::GESTURE_SUBTRACT)
		{
			RemoveBrick();
		}

		selectedBricks.anchor = selectedBricks.box;
	}

	// Reset the selected voxel aabb 
	UpdateSelectedBrick(); 
}

void WorldEditor::AddBrick()
{
	World& world = *GetWorld();
	for (int bx = selectedBricks.box.bmin[0]; bx <= selectedBricks.box.bmax[0]; bx++)
		for (int by = selectedBricks.box.bmin[1]; by <= selectedBricks.box.bmax[1]; by++)
			for (int bz = selectedBricks.box.bmin[2]; bz <= selectedBricks.box.bmax[2]; bz++)
				world.SetBrick(bx * BRICKDIM, by * BRICKDIM, bz * BRICKDIM, WHITE);
}

void WorldEditor::RemoveBrick()
{
	World& world = *GetWorld();
	for (int bx = selectedBricks.box.bmin[0]; bx <= selectedBricks.box.bmax[0]; bx++)
		for (int by = selectedBricks.box.bmin[1]; by <= selectedBricks.box.bmax[1]; by++)
			for (int bz = selectedBricks.box.bmin[2]; bz <= selectedBricks.box.bmax[2]; bz++)
				world.RemoveBrick(bx, by, bz);
}

void WorldEditor::MouseUp(int mouseButton)
{
	if ((gesture.state != GestureState::GESTURE_START && gesture.state != GestureState::GESTURE_UPDATE) ||
		gesture.button != (GestureButton)mouseButton)
	{
		return;
	}

	gesture.state = GestureState::GESTURE_POSSIBLE;
	gesture.button = GestureButton::GESTURE_DEFAULT;

	// Reset the selected voxel aabb 
	selectedBricks.box = aabb{};
	UpdateSelectedBrick();
}

// Linear to SRGB, also via https://www.shadertoy.com/view/3sfBWs
bool LessThan(float2 a, float2 b)
{
	return
		(a.x < b.x) &&
		(a.y < b.y);
}

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
	// trace the ray

	Intersection intersection; 

	if (gesture.state == GestureState::GESTURE_UPDATE || gesture.state == GestureState::GESTURE_START)
	{
		intersection = Trace(ray, tempBricks, tempGrid);
		printf("Temp Trace\n");
	}
	else 
	{
		intersection = Trace(ray);
		printf("Normal Trace\n");

	}

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
		printf("Hit World Grid - Brick X: %i   Brick Y: %i    Brick Z: %i \n", (int)brickPos.x, (int)brickPos.y, (int)brickPos.z);

	}
	else
	{

		float t = intersection.GetDistance();
		float3 N = intersection.GetNormal();
		float3 hitPoint = ray.O + ray.D * t;
		// Get position inside of the voxel to determine brick location
		float3 voxelPos = hitPoint - 0.1 * N;
		float3 brickPos = make_float3((int)voxelPos.x / BRICKDIM, (int)voxelPos.y / BRICKDIM, (int)voxelPos.z / BRICKDIM);
		
		// Adding
		if (gesture.mode == GestureMode::GESTURE_ADD)
		{
			selectedBricks.box.Grow(brickPos + N);

		}
		else if (gesture.mode == GestureMode::GESTURE_SUBTRACT)
		{
			selectedBricks.box.Grow(brickPos);
		}
		printf("Hit Brick \n");
	}

	params.selectedMin = selectedBricks.box.bmin3 * BRICKDIM;
	params.selectedMax = selectedBricks.box.bmax3 * BRICKDIM + make_float3(BRICKDIM, BRICKDIM, BRICKDIM);
}