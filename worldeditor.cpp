#include "precomp.h"
#include "worldeditor.h"


WorldEditor::WorldEditor()
{
	tempBricks = (PAYLOAD*)_aligned_malloc(CHUNKCOUNT * CHUNKSIZE, 64);
	tempGrid = (uint*)_aligned_malloc(GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * 4, 64);
	memset(tempGrid, 0, GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));
}

void WorldEditor::MouseMove(int x, int y)
{
	mousePos.x = x, mousePos.y = y;

	// just hovering
	if (gesture.state == GESTURE_POSSIBLE)
	{
		UpdateSelectedBrick();
		return;
	}

	if (gesture.state == GESTURE_START)
	{
		gesture.state = GESTURE_UPDATE;
	}

	if (gesture.button == GESTURE_LMB)
	{
		AddBrick();
	}

	if (gesture.button == GESTURE_RMB)
	{
		RemoveBrick();
	}

	UpdateSelectedBrick();
}

void WorldEditor::MouseDown(int button)
{
	if (gesture.state != GESTURE_POSSIBLE)
	{
		return;
	}

	gesture.state = GESTURE_START;
	gesture.button = (GestureButton)button;
	World& world = *GetWorld();
	memcpy(world.GetBrick(), tempBricks, BRICKSIZE * PAYLOADSIZE);
	memcpy(world.GetGrid(), tempGrid, GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * 4);

	// Delete all voxels within the selected bounding box
	if (button == GLFW_MOUSE_BUTTON_RIGHT)
	{
		RemoveBrick();
	}

	// For now, add a brick to the face of the selected voxel the cursor is on
	// Needs to be reworked when allowing for selection of multiple voxels
	if (button == GLFW_MOUSE_BUTTON_LEFT)
	{
		AddBrick();
	}

	UpdateSelectedBrick();
}

void WorldEditor::AddBrick()
{
	World& world = *GetWorld();
	float3 newBrickPos = selectedBricks.box.bmin3 + selectedBricks.N;
	world.SetBrick(newBrickPos.x * BRICKDIM, newBrickPos.y * BRICKDIM, newBrickPos.z * BRICKDIM, WHITE);
}

void WorldEditor::RemoveBrick()
{
	World& world = *GetWorld();
	for (int x = selectedBricks.box.bmin[0]; x <= selectedBricks.box.bmax[0]; x++)
		for (int y = selectedBricks.box.bmin[1]; y <= selectedBricks.box.bmax[1]; y++)
			for (int z = selectedBricks.box.bmin[2]; z <= selectedBricks.box.bmax[2]; z++)
			{
				world.SetBrick(x, y, z, 0);
			}
}

void WorldEditor::MouseUp(int button)
{
	if ((gesture.state != GESTURE_START && gesture.state != GESTURE_UPDATE) || 
		gesture.button != button)
	{
		return;
	}

	World& world = *GetWorld();
	memcpy(tempBricks, world.GetBrick(),BRICKSIZE * PAYLOADSIZE);
	memcpy(tempGrid, world.GetGrid(), GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * 4);

	gesture.state = GESTURE_POSSIBLE;
	UpdateSelectedBrick();
}

// Linear to SRGB, also via https://www.shadertoy.com/view/3sfBWs
bool LessThan(const float2 a, const float2 b)
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

	if (gesture.state == GESTURE_UPDATE || gesture.state == GESTURE_START)
	{
		intersection = Trace(ray, tempBricks, tempGrid);
		//printf("Dragging Mouse \n");
	}
	else 
	{
		intersection = Trace(ray);
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

		// Reset the selected voxel aabb 
		selectedBricks.box = aabb{};
		float3 brickPos = make_float3((int)voxelPos.x / BRICKDIM, (int)voxelPos.y / BRICKDIM, (int)voxelPos.z / BRICKDIM);
		selectedBricks.box.Grow(brickPos);
		selectedBricks.N = float3(0, 0, 0);
		printf("Hit World Grid - Brick X: %i   Brick Y: %i    Brick Z: %i \n", (int)brickPos.x, (int)brickPos.y, (int)brickPos.z);

	}
	else
	{

		float t = intersection.GetDistance();
		float3 N = intersection.GetNormal();
		float3 hitPoint = ray.O + ray.D * t;
		// Get position inside of the voxel to determine brick location
		float3 voxelPos = hitPoint - 0.1 * N;

		// Reset the selected voxel aabb 
		selectedBricks.box = aabb{};
		selectedBricks.box.Grow(make_float3((int)voxelPos.x / BRICKDIM, (int)voxelPos.y / BRICKDIM, (int)voxelPos.z / BRICKDIM));
		selectedBricks.N = N;
		printf("Hit Brick \n");
	}
}