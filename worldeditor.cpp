#include "precomp.h"
#include "worldeditor.h"

#include "lib/stb_image_write.h"

#define MAX_ALLOCATION 512000000 // max bytes allowed for undo/redo states

#define BIG_TILE_IMAGE_SCALE 8
#define TILE_IMAGE_SCALE 32

#define TILE_IMAGE_DIM BRICKDIM * TILE_IMAGE_SCALE

#pragma region NBTHelper
// Helper functions to write and read NBT tags from stream
void NBTHelper::WriteTag(std::ofstream& wf, Tag& tag)
{
	WriteTagType(wf, tag.type);

	if (tag.type == TAG_End) return;

	WriteTagName(wf, tag.name);
	if (tag.type == TAG_List)
	{
		WriteTagList(wf, tag);
	}
	if (tag.type == TAG_Compound)
	{
		WriteTagCompound(wf, tag);
	}
	if (tag.type == TAG_Int)
	{
		WriteTagInt(wf, tag);
	}
	if (tag.type == TAG_Byte_Array)
	{
		WriteTagByteArray(wf, tag);
	}
	if (tag.type == TAG_End)
	{
		int end = TAG_End;
		wf.write((char*)&end, 1);
	}
}

void NBTHelper::WriteTagList(std::ofstream& wf, Tag& tag)
{
	WriteTagType(wf, (int)tag.payload[0]);  // Write the type of tag in the list
	int numElems = tag.tags.size();
	wf.write((char*)&numElems, sizeof(int)); // Write the number of elements in the list
	for (auto tag : tag.tags)
	{
		WriteTag(wf, tag);
	}
}

void NBTHelper::WriteTagCompound(std::ofstream& wf, Tag& tag)
{
	for (auto tag : tag.tags)
	{
		WriteTag(wf, tag);
	}
}

void NBTHelper::WriteTagType(std::ofstream& wf, int tagType)
{
	wf.write((char*)&tagType, 1);
}

void NBTHelper::WriteTagInt(std::ofstream& wf, Tag& tag)
{
	wf.write((char*)&tag.payload[0], sizeof(int));
}

void NBTHelper::WriteTagName(std::ofstream& wf, std::string value)
{
	int stringLength = value.size();
	if (stringLength == 0) return;

	wf.write((char*)&stringLength, sizeof(short));
	wf.write((char*)&value[0], stringLength);
}

void NBTHelper::WriteTagByteArray(std::ofstream& wf, Tag& tag)
{
	int numBytes = tag.payload.size();
	wf.write((char*)&numBytes, sizeof(int));
	wf.write((char*)&tag.payload[0], numBytes);
}

void NBTHelper::ReadTag(std::ifstream& rf, Tag& tag, bool readTagName)
{
	ReadTagType(rf, tag.type);

	if (tag.type == TAG_End)
	{
		return;
	}

	if (readTagName) ReadTagName(rf, tag.name);

	if (tag.type == TAG_List)
	{
		ReadTagList(rf, tag);
	}
	if (tag.type == TAG_Compound)
	{
		ReadTagCompound(rf, tag);
	}
	if (tag.type == TAG_Int)
	{
		ReadTagInt(rf, tag);
	}
	if (tag.type == TAG_Byte_Array)
	{
		ReadTagByteArray(rf, tag);
	}
}

void NBTHelper::ReadTagName(std::ifstream& rf, string& tagName)
{
	int nameLength = 0;
	rf.read((char*)&nameLength, sizeof(short));
	tagName.resize(nameLength);
	rf.read((char*)&tagName[0], nameLength);
}

void NBTHelper::ReadTagType(std::ifstream& rf, int& tagType)
{
	rf.read((char*)&tagType, 1);
}

void NBTHelper::ReadTagInt(std::ifstream& rf, Tag& tag) 
{
	tag.payload.resize(sizeof(int));
	rf.read((char*)&tag.payload[0], sizeof(int));
}

void NBTHelper::ReadTagList(std::ifstream& rf, Tag& tag) 
{
	int nestedTagType = 0;
	ReadTagType(rf, nestedTagType);
	tag.payload.push_back((byte)nestedTagType);

	int numElems = 0;
	rf.read((char*)&numElems, sizeof(int)); // Write the number of elements in the list

	for (int i = 0; i < numElems; i++)
	{
		NBTHelper::Tag nestedTag;
		ReadTag(rf, nestedTag, false); // List contains unnamed tags so don't try to read name
		tag.tags.push_back(nestedTag);
	}


}

void NBTHelper::ReadTagByteArray(std::ifstream& rf, Tag& tag) 
{
	int numBytes = 0;
	rf.read((char*)&numBytes, sizeof(int));
	tag.payload.resize(numBytes);
	rf.read((char*)&tag.payload[0], numBytes);
};

void NBTHelper::ReadTagCompound(std::ifstream& rf, Tag& tag) 
{
	while (true)
	{
		NBTHelper::Tag nestedTag;
		ReadTag(rf, nestedTag);
		tag.tags.push_back(nestedTag);
		if (nestedTag.type == TAG_End)
			return;
	}
};
#pragma endregion



#pragma region Initialization
// Simple helper function to load an image into a OpenGL texture with common settings
bool LoadTextureFromFile(const char* filename, GLuint* out_texture, int* out_width, int* out_height)
{
	// Load from file
	int image_width = 0;
	int image_height = 0;
	unsigned char* image_data = stbi_load(filename, &image_width, &image_height, NULL, 4);
	if (image_data == NULL)
		return false;

	// Create a OpenGL texture identifier
	GLuint image_texture;
	glGenTextures(1, &image_texture);
	glBindTexture(GL_TEXTURE_2D, image_texture);

	// Setup filtering parameters for display
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	// Upload pixels into texture
#if defined(GL_UNPACK_ROW_LENGTH) && !defined(__EMSCRIPTEN__)
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
#endif
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image_width, image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);
	stbi_image_free(image_data);

	*out_texture = image_texture;
	*out_width = image_width;
	*out_height = image_height;

	return true;
}

float blueNoiseSampler(const uint* blueNoise, int x, int y, int sampleIndex, int sampleDimension)
{
	// Adapated from E. Heitz. Arguments:
	// sampleIndex: 0..255
	// sampleDimension: 0..255
	x &= 127, y &= 127, sampleIndex &= 255, sampleDimension &= 255;
	// xor index based on optimized ranking
	int rankedSampleIndex = (sampleIndex ^ blueNoise[sampleDimension + (x + y * 128) * 8 + 65536 * 3]) & 255;
	// fetch value in sequence
	int value = blueNoise[sampleDimension + rankedSampleIndex * 256];
	// if the dimension is optimized, xor sequence value based on optimized scrambling
	value ^= blueNoise[(sampleDimension & 7) + (x + y * 128) * 8 + 65536];
	// convert to float and return
	float retVal = (0.5f + value) * (1.0f / 256.0f) /* + noiseShift (see LH2) */;
	if (retVal >= 1) retVal -= 1;
	return retVal;
}

void WorldEditor::LoadAssets()
{
	const std::filesystem::path tiles{ "assets/Tiles" };
	const std::filesystem::path bigTiles{ "assets/BigTiles" };
	const std::filesystem::path sprites{ "assets/Sprites" };

	TileManager* tileManager = TileManager::GetTileManager();
	World& world = *GetWorld();
	RenderParams& params = world.GetRenderParams();

	const uchar* data8 = (const uchar*)sob256_64; // tables are 8 bit per entry
	uint* blueNoise = new uint[65536 * 5]; // we want a full uint per entry
	for (int i = 0; i < 65536; i++) blueNoise[i] = data8[i]; // convert
	data8 = (uchar*)scr256_64;
	for (int i = 0; i < (128 * 128 * 8); i++) blueNoise[i + 65536] = data8[i];
	data8 = (uchar*)rnk256_64;
	for (int i = 0; i < (128 * 128 * 8); i++) blueNoise[i + 3 * 65536] = data8[i];


	auto CreatePreview = [&blueNoise, &world](const char* filename, float3 O, float3 D)
	{
		Surface* surface = new Surface(TILE_IMAGE_DIM, TILE_IMAGE_DIM);

		mat4 M = mat4::LookAt(O, O + D);
		float3 cam = TransformPosition(make_float3(0), M);
		float3 p0 = TransformPosition(make_float3(-1, 1, 2), M);
		float3 p1 = TransformPosition(make_float3(1, 1, 2), M);
		float3 p2 = TransformPosition(make_float3(-1, -1, 2), M);

		int numSamplesPerPixel = 5;
		for (int y = 0; y < TILE_IMAGE_DIM; y++)
		{
			for (int x = 0; x < TILE_IMAGE_DIM; x++)
			{
				float3 color = make_float3(0.0f, 0.0f, 0.0f);

				uint val = 0;
				for (int s = 0; s < numSamplesPerPixel; s++)
				{
					float u = (float(x) + RandomFloat()) / static_cast<float>(TILE_IMAGE_DIM);
					float v = (float(y) + RandomFloat()) / static_cast<float>(TILE_IMAGE_DIM);

					Ray ray;
					ray.O = O;
					ray.D = normalize((p0 + (p1 - p0) * u + (p2 - p0) * v) - O);
					ray.t = 1e34f;

					Intersection intersection = Trace(ray);
					if (intersection.GetVoxel())
					{
						color = ToFloatRGB(intersection.GetVoxel());
						const float3 brdf = color * INVPI;

						float3 incoming = make_float3(0.0f, 0.0f, 0.0f);
						const float skyLightScale = 3.0f;
						const float r0 = blueNoiseSampler(blueNoise, x, y, 0, 0);
						const float r1 = blueNoiseSampler(blueNoise, x, y, 0, 1);
						const float3 R = DiffuseReflectionCosWeighted(r0, r1, intersection.GetNormal());
						uint side2;
						float dist2;
						const float3 shadingPoint = D * intersection.GetDistance() + ray.O;
						const float3 shadingPointOffset = shadingPoint + 0.1 * intersection.GetNormal();
						Ray ray2;
						ray2.O = shadingPointOffset;
						ray2.D = R;
						ray2.t = 1e34f;

						Intersection intersection2 = Trace(ray2);
						const float3 N2 = intersection2.GetNormal();
						float3 toAdd = make_float3(skyLightScale), M = intersection.GetNormal();
						// TO-DO: Add alpha blend here
						if (intersection2.GetVoxel() != 0)
							toAdd *= INVPI * ToFloatRGB(intersection2.GetVoxel()), M = N2;

						float4* skyLight = world.GetSkyLight();
						float4 sky;
						if (M.x < -0.9f) sky = skyLight[0];
						if (M.x > 0.9f) sky = skyLight[1];
						if (M.y < -0.9f) sky = skyLight[2];
						if (M.y > 0.9f) sky = skyLight[3];
						if (M.z < -0.9f) sky = skyLight[4];
						if (M.z > 0.9f) sky = skyLight[5];
						incoming += toAdd * make_float3(sky.x, sky.y, sky.z);

						color += incoming * brdf;
						val |= 255 << 24;
					}
				}

				float scale = 1.0f / numSamplesPerPixel;
				val |= ((std::min<uint32_t>(static_cast<uint32_t>(sqrt(scale * color.x) * 255), 255) << 0) |
					(std::min<uint32_t>(static_cast<uint32_t>(sqrt(scale * color.y) * 255), 255) << 8) |
					(std::min<uint32_t>(static_cast<uint32_t>(sqrt(scale * color.z) * 255), 255) << 16));

				surface->Plot((TILE_IMAGE_DIM - 1 - x), y, val);
			}
		}

		stbi_write_png(filename, TILE_IMAGE_DIM, TILE_IMAGE_DIM, 4, surface->buffer, TILE_IMAGE_DIM * sizeof(uint32_t));
	};

	float3 cameraOrigin = make_float3(12, 10, -8);
	// Load any saved Tiles
	for (auto const& dir_entry : std::filesystem::directory_iterator{ tiles })
	{
		std::filesystem::path p = dir_entry.path();
		if (p.extension() == ".vox")
		{
			// Load the tile
			uint tileIdx = LoadTile(dir_entry.path().string().c_str());

			// Check to see if the tile has a preview png already
			std::filesystem::path png = ".png";
			p.replace_extension(png);
			if (!std::filesystem::exists(p)) // Construct a preview for the tile
			{
				PAYLOAD* voxels = tileManager->tile[tileIdx]->voxels;

				world.DrawTile(tileIdx, 0, 0, 0);

				float3 O = cameraOrigin;
				float3 D = make_float3(-0.5, -0.5, 0.7);
				CreatePreview(p.string().c_str(), O, D);
			}

			int my_image_width = 0;
			int my_image_height = 0;
			GLuint my_image_texture = 0;

			bool ret = LoadTextureFromFile(p.string().c_str(), &my_image_texture, &my_image_width, &my_image_height);
			IM_ASSERT(ret);
			loadedTiles.push_back(std::make_pair(tileIdx, my_image_texture));
		}
	}

	// Load any saved BigTiles
	for (auto const& dir_entry : std::filesystem::directory_iterator{ bigTiles })
	{
		std::filesystem::path p = dir_entry.path();
		if (p.extension() == ".vox")
		{
			// Load the tile
			uint bigTileIdx = LoadBigTile(dir_entry.path().string().c_str());

			// Check to see if the tile has a preview png already
			std::filesystem::path png = ".png";
			p.replace_extension(png);
			if (!std::filesystem::exists(p)) // Construct a preview for the tile
			{
				
				world.DrawBigTile(bigTileIdx, 0, 0, 0);

				float3 O = cameraOrigin * 2;
				float3 D = make_float3(-0.5, -0.5, 0.7);
				CreatePreview(p.string().c_str(), O, D);
			}
			int my_image_width = 0;
			int my_image_height = 0;
			GLuint my_image_texture = 0;

			bool ret = LoadTextureFromFile(p.string().c_str(), &my_image_texture, &my_image_width, &my_image_height);
			IM_ASSERT(ret);
			loadedBigTiles.push_back(std::make_pair(bigTileIdx, my_image_texture));
		}
	}

	// Load any saved BigTiles
	for (auto const& dir_entry : std::filesystem::directory_iterator{ sprites })
	{
		std::filesystem::path p = dir_entry.path();
		if (p.extension() == ".vox")
		{
			// Load the tile
			uint spriteIdx = LoadSprite(dir_entry.path().string().c_str());

			// Check to see if the tile has a preview png already
			std::filesystem::path png = ".png";
			p.replace_extension(png);
			if (!std::filesystem::exists(p)) // Construct a preview for the tile
			{

				world.StampSpriteTo(spriteIdx, 0, 0, 0);
				Sprite* sprite = world.GetSpriteList()[spriteIdx];
				int3 spriteSize = sprite->frame[0]->size;
				float3 cameraShift = make_float3(spriteSize) / 8.0f;
				cameraShift.x = max(cameraShift.x, 1.0f), cameraShift.y = max(cameraShift.y, 1.0f), cameraShift.z = max(cameraShift.x, 1.0f);
				float3 O = cameraOrigin * cameraShift;
				float3 D = make_float3(-0.5, -0.5, 0.7);
				CreatePreview(p.string().c_str(), O, D);

				world.Clear();
			}
			int my_image_width = 0;
			int my_image_height = 0;
			GLuint my_image_texture = 0;

			bool ret = LoadTextureFromFile(p.string().c_str(), &my_image_texture, &my_image_width, &my_image_height);
			IM_ASSERT(ret);
			loadedSprites.push_back(std::make_pair(spriteIdx, my_image_texture));
		}
	}

	delete[] blueNoise;
	world.Clear();
}

WorldEditor::WorldEditor()
{
	tempBricks = (PAYLOAD*)_aligned_malloc(CHUNKCOUNT * CHUNKSIZE, 64);
	tempGrid = (uint*)_aligned_malloc(GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * 4, 64);
	tempZeroes = (uint*)_aligned_malloc(GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint), 64);
	stateHead = (State*)calloc(1, sizeof(State));
	stateCurrent = stateTail = stateHead;

	memset(tempBricks, 0, CHUNKCOUNT * CHUNKSIZE);
	memset(tempGrid, 0, GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));
	memset(tempZeroes, BRICKSIZE, GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));

	LoadAssets();
}
#pragma endregion 



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
		selected.box = aabb{};
		UpdateSelectedBox();
		return;
	}

	// Don't allow multi add/remove for sprites
	if (gesture.mode & GestureMode::GESTURE_MULTI && gesture.size != GestureSize::GESTURE_SPRITE)
	{
		MultiAddRemove();
	}
	else 
	{
		auto oldBox = selected.box;
		selected.box = aabb{};
		UpdateSelectedBox();

		// Ignore when mouse hovers over the same brick to avoid multiple add/removes
		if (oldBox == selected.box) return;

		float3 boxPos = selected.box.bmin3;
		if (gesture.mode & GestureMode::GESTURE_REMOVE) Remove(boxPos.x, boxPos.y, boxPos.z);
		else Add(boxPos.x, boxPos.y, boxPos.z);
	}
	UpdateSelectedBox();
}

void WorldEditor::KeyDown(int key)
{
	if (!enabled) return;

	switch (key) {
	case GLFW_KEY_LEFT_CONTROL:
		selectedKeys |= GestureKey::GESTURE_CTRL;
		UpdateGestureMode();
		selected.box = aabb{};
		UpdateSelectedBox();
		break;
	case GLFW_KEY_LEFT_SHIFT:
		selectedKeys |= GestureKey::GESTURE_SHIFT;
		UpdateGestureMode();
		selected.box = aabb{};
		UpdateSelectedBox();
		break;
	case GLFW_KEY_Z:
		if (selectedKeys == GestureKey::GESTURE_CTRL)
		{
			Undo();
			selected.box = aabb{};
			UpdateSelectedBox();
		}
		break;
	case GLFW_KEY_Y:
		if (selectedKeys == GestureKey::GESTURE_CTRL)
		{
			Redo();
			selected.box = aabb{};
			UpdateSelectedBox();
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
		selected.box = aabb{};
		UpdateSelectedBox();
		break;
	case GLFW_KEY_LEFT_SHIFT:
		selectedKeys ^= GestureKey::GESTURE_SHIFT;
		UpdateGestureMode();
		selected.box = aabb{};
		UpdateSelectedBox();
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
		editedBricks.clear();
		gesture.state = GestureState::GESTURE_ACTIVE;

		// Save a backup version of the grid/brick/info 
		World& world = *GetWorld();
		memcpy(tempBricks, world.GetBrick(), CHUNKCOUNT * CHUNKSIZE);
		memcpy(tempGrid, world.GetGrid(), GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));
		memcpy(tempZeroes, world.GetZeroes(), GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));

		float3 brickPos = selected.box.bmin3;
		// Add or remove the intial brick where the mouse is hovered
		if (gesture.mode & GestureMode::GESTURE_REMOVE) Remove(brickPos.x, brickPos.y, brickPos.z);
		else Add(brickPos.x, brickPos.y, brickPos.z);

		// Set our anchor for multi add/remove
		selected.anchor = selected.box;
	}

	UpdateSelectedBox();
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
		selected.box = aabb{};
		UpdateSelectedBox();
	}
}
#pragma endregion



#pragma region Editing
// Update what bricks have been edited during this gesture depending on the gesture type
void WorldEditor::UpdateEditedBricks(uint x, uint y, uint z)
{
	if (gesture.size == GestureSize::GESTURE_BIG_TILE)
	{
		for (uint _x = 0; _x < 2; _x++)
			for (uint _y = 0; _y < 2; _y++)
				for (uint _z = 0; _z < 2; _z++)
					editedBricks.insert(make_int3(x / BRICKDIM + _x, y / BRICKDIM + _y, z / BRICKDIM + _z));
	}
	else if (gesture.size == GestureSize::GESTURE_SPRITE)
	{
		World& world = *GetWorld();
		int3 spriteSize = world.GetSpriteList()[selectedSpriteIdx]->frame[0]->size;
		for (uint _x = 0; _x < spriteSize.x; _x++)
			for (uint _y = 0; _y < spriteSize.y; _y++)
				for (uint _z = 0; _z < spriteSize.z; _z++)
					editedBricks.insert(make_int3((x + _x) / BRICKDIM, (y + _y) / BRICKDIM, (z + _z) / BRICKDIM));
	}
	else
	{
		editedBricks.insert((make_int3)(x / BRICKDIM, y / BRICKDIM, z / BRICKDIM));
	}

	return;
}

// Get the scale required to draw the selected box
int3 WorldEditor::GetBoxScale()
{
	switch (gesture.size)
	{
	case GestureSize::GESTURE_TILE:
	case GestureSize::GESTURE_BRICK:
		return make_int3(BRICKDIM);
	case GestureSize::GESTURE_BIG_TILE:
		return make_int3(BRICKDIM * 2);
	case GestureSize::GESTURE_SPRITE:
	{
		World& world = *GetWorld();
		int3 spriteSize = world.GetSpriteList()[selectedSpriteIdx]->frame[0]->size;
		return spriteSize;
	}
	case GestureSize::GESTURE_VOXEL:
	default:
		return make_int3(1);
	}
}

// Add back lights in the world after a Redo/Undo
void WorldEditor::AddBackLights(uint bx, uint by, uint bz)
{
	World& world = *GetWorld();
	for (uint ly = 0; ly < BRICKDIM; ly++)
		for (uint lz = 0; lz < BRICKDIM; lz++)
			for (uint lx = 0; lx < BRICKDIM; lx++)
			{
				uint v = world.Get(lx + bx * BRICKDIM, ly + by * BRICKDIM, lz + bz * BRICKDIM);

				if (IsEmitter(v))
				{
					world.AddLight(int3(lx + bx * BRICKDIM, ly + by * BRICKDIM, lz + bz * BRICKDIM), 1, v);
				}
			}
}

// Allow the adding/removing of multiple bricks in the seclected box
void WorldEditor::MultiAddRemove()
{
	aabb oldBox = selected.box;
	selected.box = selected.anchor;
	UpdateSelectedBox();

	if (oldBox == selected.box) return;

	World& world = *GetWorld();
	uint* brick = world.GetBrick();
	uint* grid = world.GetGrid();

	aabb newBox = selected.box;
	// Compute the overlap betweeen the old and new aabbs

	int3 boxScale = GetBoxScale();

	oldBox.bmin3 = make_float3(make_int3(oldBox.bmin3) / boxScale);
	oldBox.bmax3 = make_float3(make_int3(oldBox.bmax3) / boxScale);
	newBox.bmin3 = make_float3(make_int3(newBox.bmin3) / boxScale);
	newBox.bmax3 = make_float3(make_int3(newBox.bmax3) / boxScale);
	aabb intersection = oldBox.Intersection(newBox);

	auto RestoreTempBrickVal = [&](uint bx, uint by, uint bz) {
		const uint cellIndex = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
		const uint tempCellVal = tempGrid[cellIndex];
		const uint curCellVal = grid[cellIndex];

		// If we're reintroducing what was an empty/solid brick, remove the current one
		if (IsSolidGridCell(tempCellVal))
		{
			world.RemoveBrick(bx, by, bz);
			world.AddBrick(bx, by, bz, tempCellVal >> 1);
			return;
		}

		// Get current and new brick offsets
		uint tempBrickBufferOffset = (tempCellVal >> 1);
		uint curBrickBufferOffset = (curCellVal >> 1);

		uint brickBufferOffset = tempBrickBufferOffset;
		uint cellVal = tempCellVal;

		// If the current brick is empty/solid we need to create a NewBrick
		if (IsSolidGridCell(curCellVal)) brickBufferOffset = world.NewBrick(), cellVal = (brickBufferOffset << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(brick + brickBufferOffset * BRICKSIZE, tempBricks + tempBrickBufferOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
		grid[cellIndex] = cellVal;
		world.GetZeroes()[brickBufferOffset] = tempZeroes[tempBrickBufferOffset];
		world.Mark(brickBufferOffset);
		AddBackLights(bx, by, bz);
	};

	if (gesture.size == GestureSize::GESTURE_BRICK || gesture.size == GestureSize::GESTURE_TILE)
	{
		// We've updated the selected box so now it doesn't cover the same area
		// For all individual updated values, put back the old ones
		for (int x = oldBox.bmin3.x; x <= oldBox.bmax3.x; x++)
			for (int y = oldBox.bmin3.y; y <= oldBox.bmax3.y; y++)
				for (int z = oldBox.bmin3.z; z <= oldBox.bmax3.z; z++)
				{
					__m128 b4 = _mm_setr_ps(x, y, z, 0);
					if (intersection.Contains(b4))
						continue;

					RestoreTempBrickVal(x, y, z);
				}
	}

	if (gesture.size == GESTURE_BIG_TILE)
	{
		// We've updated the selected box so now it doesn't cover the same area
		// For all individual updated values, put back the old ones
		for (int x = oldBox.bmin3.x; x <= oldBox.bmax3.x; x++)
			for (int y = oldBox.bmin3.y; y <= oldBox.bmax3.y; y++)
				for (int z = oldBox.bmin3.z; z <= oldBox.bmax3.z; z++)
				{
					__m128 b4 = _mm_setr_ps(x, y, z, 0);
					if (intersection.Contains(b4))
						continue;

					for (int _x = 0; _x < 2; _x++)
						for (int _y = 0; _y < 2; _y++)
							for (int _z = 0; _z < 2; _z++)
								RestoreTempBrickVal(x * 2 + _x, y * 2 + _y, z * 2 + _z);
				}
	}

	if (gesture.size == GESTURE_VOXEL)
	{
		// We've updated the selected box so now it doesn't cover the same area
		// For all individual updated voxels, put back the old ones
		for (int x = oldBox.bmin3.x; x <= oldBox.bmax3.x; x++)
			for (int y = oldBox.bmin3.y; y <= oldBox.bmax3.y; y++)
				for (int z = oldBox.bmin3.z; z <= oldBox.bmax3.z; z++)
				{
					__m128 b4 = _mm_setr_ps(x, y, z, 0);
					if (intersection.Contains(b4))
						continue;

					const uint bx = x / BRICKDIM;
					const uint by = y / BRICKDIM;
					const uint bz = z / BRICKDIM;
					if (bx >= GRIDWIDTH || by >= GRIDHEIGHT || bz >= GRIDDEPTH) return;
					const uint cellIndex = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
					// obtain current brick identifier from top-level grid
					uint cellValue = grid[cellIndex], brickBufferOffset = cellValue >> 1;

					// calculate the position of the voxel inside the brick
					const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
					const uint voxelIndex = brickBufferOffset * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM;
					const uint oldVoxelValue = tempBricks[voxelIndex];
					world.Set(x, y, z, oldVoxelValue);
				}
	}


	// Clear all edited bricks and fill it with bricks form the new selected box aabb
	editedBricks.clear();
	// Either add or remove the values in the new selected box aabb
	for (int x = newBox.bmin3.x; x <= newBox.bmax3.x; x++)
		for (int y = newBox.bmin3.y; y <= newBox.bmax3.y; y++)
			for (int z = newBox.bmin3.z; z <= newBox.bmax3.z; z++)
			{
				__m128 b4 = _mm_setr_ps(x, y, z, 0);
				if (intersection.Contains(b4))
				{
					UpdateEditedBricks(x * boxScale.x, y * boxScale.y, z * boxScale.z);
					continue;
				}

				if (gesture.mode & GestureMode::GESTURE_REMOVE)
				{
					Remove(x * boxScale.x, y * boxScale.y, z * boxScale.z);
				}
				else
				{
					Add(x * boxScale.x, y * boxScale.y, z * boxScale.z);
				}
			}
}

// Add depending on the selected gesture
void WorldEditor::Add(uint vx, uint vy, uint vz)
{
	World& world = *GetWorld();

	// Convert voxel xyz to brick coords
	uint bx = vx / BRICKDIM;
	uint by = vy / BRICKDIM;
	uint bz = vz / BRICKDIM;

	if (gesture.size == GestureSize::GESTURE_VOXEL)
	{
		world.Set(vx, vy, vz, voxelValue);
		UpdateEditedBricks(vx, vy, vz);
	}
	else 
	{
		// Avoid adding a brick already added during this gesture
		int3 brickPos = make_int3(bx, by, bz);
		if (!editedBricks.count(brickPos))
		{
			UpdateEditedBricks(vx, vy, vz);
			if (gesture.size == GestureSize::GESTURE_TILE)
			{
				world.DrawTile(selectedTileIdx, bx, by, bz);
			}
			else if (gesture.size == GestureSize::GESTURE_BIG_TILE)
			{
				world.DrawBigTile(selectedBigTileIdx, bx / 2, by / 2, bz / 2);
			}
			else if (gesture.size == GestureSize::GESTURE_BRICK)
			{
				world.AddBrick(bx, by, bz, voxelValue);
			}
			else if (gesture.size == GestureSize::GESTURE_SPRITE)
			{
				world.StampSpriteTo(selectedSpriteIdx, vx, vy, vz);
			}
		}
	}
}

// Remove depending on the selected gesture
void WorldEditor::Remove(uint vx, uint vy, uint vz)
{
	World& world = *GetWorld();

	uint bx = vx / BRICKDIM;
	uint by = vy / BRICKDIM;
	uint bz = vz / BRICKDIM;

	if (gesture.size == GestureSize::GESTURE_VOXEL)
	{
		world.Set(vx, vy, vz, 0);
		UpdateEditedBricks(vx, vy, vz);
	}
	else
	{
		// Avoid removing a brick already removed during this gesture
		int3 brickPos = make_int3(bx, by, bz);
		if (!editedBricks.count(brickPos))
		{
			UpdateEditedBricks(vx, vy, vz);
			if (gesture.size == GestureSize::GESTURE_TILE || gesture.size == GestureSize::GESTURE_BRICK)
			{
				world.RemoveBrick(bx, by, bz);
			}
			else if (gesture.size == GestureSize::GESTURE_BIG_TILE)
			{
				world.RemoveBigTile(bx, by, bz);
			}
			else if (gesture.size == GestureSize::GESTURE_SPRITE)
			{
				int3 spriteSize = world.GetSpriteList()[selectedSpriteIdx]->frame[0]->size;
				for (uint _x = 0; _x < spriteSize.x; _x++)
					for (uint _y = 0; _y < spriteSize.y; _y++)
						for (uint _z = 0; _z < spriteSize.z; _z++)
							world.Set(vx + _x, vy + _y, vz + _z, 0);
			}
		}
	}

}

// Update which brick is currently selected by the mouse cursor 
void WorldEditor::UpdateSelectedBox()
{
	World& world = *GetWorld();
	RenderParams& params = world.GetRenderParams();

	int3 boxScale = GetBoxScale();

	const float2 uv = make_float2(mousePos.x * params.oneOverRes.x, mousePos.y * params.oneOverRes.y);
	const float3 P = params.p0 + (params.p1 - params.p0) * uv.x + (params.p2 - params.p0) * uv.y;

	Ray ray;
	ray.O = params.E;
	ray.D = normalize(P - params.E);
	ray.t = 1e34f;

	// Trace the ray using the previous grid/brick state if gesture is active
	Intersection intersection = (gesture.state == GestureState::GESTURE_ACTIVE) ? Trace(ray, tempBricks, tempGrid) : Trace(ray);

	// Check to see if we hit the world grid using method from  Majercik et al. https://jcgt.org/published/0007/03/04/paper-lowres.pdf
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


		float3 hitPoint = ray.O + ray.D * distance;
		// Get position inside of the grid to determine brick location
		float3 gridPos = hitPoint + 0.5 * normal;

		if (gridPos.x < 0 || gridPos.y < 0 || gridPos.z < 0 || gridPos.x > MAPWIDTH || gridPos.y > MAPHEIGHT || gridPos.z > MAPDEPTH)
			return;

		float3 minBoxPos;
		if (gesture.size == GestureSize::GESTURE_SPRITE)
		{
			// Allow sprites to be placed anywhere
			minBoxPos = make_float3(make_int3(gridPos));
		}
		else
		{
			minBoxPos = make_float3(make_int3(gridPos / make_float3(boxScale)) * boxScale);
		}

		selected.box.Grow(minBoxPos);
		float3 maxBoxPos = minBoxPos + make_float3(boxScale - make_int3(1));
		selected.box.Grow(maxBoxPos);
	}
	else
	{
		float t = intersection.GetDistance();
		float3 normal = intersection.GetNormal();
		float3 hitPoint = ray.O + ray.D * t;

		// Get position inside of the voxel to determine brick location
		float3 voxelPos = hitPoint - 0.1 * normal;
		float3 minBoxPos;

		if (gesture.size == GestureSize::GESTURE_SPRITE)
		{
			// Allow sprites to be placed anywhere
			minBoxPos = make_float3(make_int3(voxelPos));
		}
		else
		{
			minBoxPos = make_float3(make_int3(voxelPos / make_float3(boxScale)) * boxScale);
		}

		if (gesture.mode & GestureMode::GESTURE_REMOVE)
		{
			selected.box.Grow(minBoxPos);
			selected.box.Grow(minBoxPos + make_float3(boxScale - make_int3(1)));

		}
		else
		{
			float3 newBoxPos;
			if (gesture.size == GestureSize::GESTURE_SPRITE)
			{
				newBoxPos = minBoxPos + normal;
			}
			else
			{
				newBoxPos = minBoxPos + normal * boxScale;

			}

			// If the normal puts us outside of the grid, just use the original intersection
			if (newBoxPos.x < 0 || newBoxPos.y < 0 || newBoxPos.z < 0 || newBoxPos.x > MAPWIDTH || newBoxPos.y > MAPHEIGHT || newBoxPos.z > MAPDEPTH)
			{
				selected.box.Grow(minBoxPos);
				selected.box.Grow(minBoxPos + make_float3(boxScale - make_int3(1)));
			}
			else
			{
				selected.box.Grow(newBoxPos);
				selected.box.Grow(newBoxPos + make_float3(boxScale - make_int3(1)));
			}
		}
	}

	// Update rendering params to trace the selected box outline
	params.selectedMin = selected.box.bmin3;
	params.selectedMax = selected.box.bmax3;
	if (gesture.size == GestureSize::GESTURE_BIG_TILE || gesture.size == GestureSize::GESTURE_SPRITE) params.wireBoxWidth = 0.5f;
	if (gesture.size == GestureSize::GESTURE_TILE || gesture.size == GestureSize::GESTURE_BRICK) params.wireBoxWidth = 0.3f;
	if (gesture.size == GestureSize::GESTURE_VOXEL) params.wireBoxWidth = 0.1f;

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
#pragma endregion



#pragma region State
// Reset the world editor so when re-enabled we're clean 
void WorldEditor::ResetEditor()
{
	gesture = Gesture{};
	selected = Selected{};
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
	uint* zeroes = world.GetZeroes();

	int numBricks = stateCurrent->numBricks;
	for (int idx = 0; idx < numBricks; idx++)
	{
		int3 b = stateCurrent->editedBricks[idx];
		const uint cellIndex = b.x + b.z * GRIDWIDTH + b.y * GRIDWIDTH * GRIDDEPTH;

		const uint oldCellValue = stateCurrent->oldCellValues[idx];
		const uint newCellValue = stateCurrent->newCellValues[idx];

		// If we're restoring what was an empty/solid brick, remove the current one
		if (IsSolidGridCell(oldCellValue))
		{
			world.RemoveBrick(b.x, b.y, b.z);
			world.AddBrick(b.x, b.y, b.z, oldCellValue >> 1);
			continue;
		}

		// Get old and new brick offsets (g1)
		uint oldbrickBufferOffset = (oldCellValue >> 1);
		uint newbrickBufferOffset = (newCellValue >> 1);

		uint brickBufferOffset = oldbrickBufferOffset;
		uint cellValue = oldCellValue;

		// If the current brick is empty/solid we need to create a NewBrick
		if (IsSolidGridCell(newCellValue)) brickBufferOffset = world.NewBrick(), cellValue = (brickBufferOffset << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(bricks + (brickBufferOffset * BRICKSIZE), stateCurrent->oldBricks + (idx * BRICKSIZE), BRICKSIZE * PAYLOADSIZE);
		grid[cellIndex] = cellValue;
		zeroes[brickBufferOffset] = stateCurrent->oldBrickZeroes[idx];
		world.Mark(brickBufferOffset);
		AddBackLights(b.x, b.y, b.z);
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
	uint* zeroes = world.GetZeroes();

	int numBricks = stateCurrent->numBricks;
	for (int idx = 0; idx < numBricks; idx++)
	{
		int3 b = stateCurrent->editedBricks[idx];
		const uint cellIndex = b.x + b.z * GRIDWIDTH + b.y * GRIDWIDTH * GRIDDEPTH;

		const uint oldCellValue = stateCurrent->oldCellValues[idx];
		const uint newCellValue = stateCurrent->newCellValues[idx];

		// If we're restoring what was an empty/solid brick, just remove the current one
		if (IsSolidGridCell(newCellValue))
		{
			world.RemoveBrick(b.x, b.y, b.z);
			world.AddBrick(b.x, b.y, b.z, newCellValue >> 1);
			continue;
		}

		// Get old and new brick offsets (g1)
		uint oldbrickBufferOffset = (oldCellValue >> 1);
		uint newbrickBufferOffset = (newCellValue >> 1);

		uint brickBufferOffset = newbrickBufferOffset;
		uint cellValue = newCellValue;

		// If the current brick is empty/solid we need to create a NewBrick
		if (IsSolidGridCell(oldCellValue)) brickBufferOffset = world.NewBrick(), cellValue = (brickBufferOffset << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(bricks + (brickBufferOffset * BRICKSIZE), stateCurrent->newBricks + (idx * BRICKSIZE), BRICKSIZE * PAYLOADSIZE);
		grid[cellIndex] = cellValue;
		zeroes[brickBufferOffset] = stateCurrent->oldBrickZeroes[idx];
		world.Mark(brickBufferOffset);
		AddBackLights(b.x, b.y, b.z);
	}
}

// Called after a gesture has been completed so we can update the history
void WorldEditor::SaveState()
{
	World& world = *GetWorld();
	uint* grid = world.GetGrid();
	PAYLOAD* bricks = world.GetBrick();
	uint* zeroes = world.GetZeroes();
	uint* trash = world.GetTrash();

	if (!CreateNewState())
	{
		undoEnabled = false;
		return;
	}

	int numBricks = editedBricks.size();
	stateCurrent->numBricks = numBricks;
	stateCurrent->newBricks = (PAYLOAD*)_aligned_malloc(numBricks * BRICKSIZE * PAYLOADSIZE, 64);
	stateCurrent->oldBricks = (PAYLOAD*)_aligned_malloc(numBricks * BRICKSIZE * PAYLOADSIZE, 64);
	stateCurrent->newCellValues = (uint*)_aligned_malloc(numBricks * 4, 64);
	stateCurrent->oldCellValues = (uint*)_aligned_malloc(numBricks * 4, 64);
	stateCurrent->newBrickZeroes = (uint*)_aligned_malloc(numBricks * sizeof(uint), 64);
	stateCurrent->oldBrickZeroes = (uint*)_aligned_malloc(numBricks * sizeof(uint), 64);
	stateCurrent->editedBricks = (int3*)_aligned_malloc(numBricks * sizeof(int3), 64);

	int index = 0;
	for (auto itr = editedBricks.begin(); itr != editedBricks.end(); itr++)
	{
		int3 brick = *itr;
		const uint cellIndex = brick.x + brick.z * GRIDWIDTH + brick.y * GRIDWIDTH * GRIDDEPTH;

		const uint tempCellValue = tempGrid[cellIndex];
		const uint currCellValue = grid[cellIndex];

		uint tempBrickBufferOffset = (tempCellValue >> 1);
		uint brickBufferOffset = (currCellValue >> 1);

		// Old brick is empty/or a solid brick
		if (IsSolidGridCell(tempCellValue))
		{
			memset(stateCurrent->oldBricks + (index * BRICKSIZE), 0, BRICKSIZE * PAYLOADSIZE);
		}
		else
		{
			memcpy(stateCurrent->oldBricks + (index * BRICKSIZE), tempBricks + tempBrickBufferOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
		}
		
		// New brick is empty/or a solid brick
		if (IsSolidGridCell(currCellValue))
		{
			memset(stateCurrent->newBricks + (index * BRICKSIZE), 0, BRICKSIZE * PAYLOADSIZE);
		}
		else
		{
			memcpy(stateCurrent->newBricks + (index * BRICKSIZE), bricks + brickBufferOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
			stateCurrent->oldBrickZeroes[index] = tempZeroes[tempBrickBufferOffset];
			stateCurrent->newBrickZeroes[index] = zeroes[brickBufferOffset];
		}

		stateCurrent->oldCellValues[index] = tempCellValue;
		stateCurrent->newCellValues[index] = currCellValue;


		stateCurrent->editedBricks[index] = brick;
		index++;
	}


	// Add the amount of allocated memory 
	allocatedUndo += 2 * (numBricks * BRICKSIZE * PAYLOADSIZE);
	allocatedUndo += 2 * (numBricks * 4);
	allocatedUndo += 2 * (numBricks * sizeof(uint));
	allocatedUndo += (numBricks * sizeof(int3));


	// Check to see if the state has actually changed....if not don't add it to the history
	if (memcmp(stateCurrent->oldCellValues, stateCurrent->newCellValues, numBricks * 4) == 0 && 
		memcmp(stateCurrent->oldBricks, stateCurrent->newBricks, numBricks * BRICKSIZE * PAYLOADSIZE) == 0)
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
	_aligned_free(state->oldCellValues);
	_aligned_free(state->newCellValues);
	_aligned_free(state->oldBrickZeroes);
	_aligned_free(state->newBrickZeroes);
	_aligned_free(state->editedBricks);
	free(state);

	allocatedUndo -= 2 * (numBricks * BRICKSIZE * PAYLOADSIZE);
	allocatedUndo -= 2 * (numBricks * 4);
	allocatedUndo -= 2 * (numBricks * sizeof(uint));
	allocatedUndo -= (numBricks * sizeof(int3));
}

void WorldEditor::LoadWorld()
{
	using namespace NBTHelper;

	World& world = *GetWorld();
	world.Clear();

	unsigned int* brick = world.GetBrick();
	uint* grid = world.GetGrid();
	uint* zeroes = world.GetZeroes();

	ifstream rf("worlddata.nbt", ios::in | ios::binary);
	if (!rf)
	{
		printf("Error opening file worlddata.nbt\n");
		return;
	}

	// Populate our tags from the worlddata.nbt file
	std::vector<Tag> topLevelTags;
	while (rf.good())
	{
		Tag tag;
		ReadTag(rf, tag);
		topLevelTags.push_back(tag);
	}
	rf.close();

	// Iterate through read tags and populate the world data
	for (auto tag : topLevelTags)
	{
		if (tag.name == "grid list")
		{
			for (auto compoundTag : tag.tags)
			{
				uint value = 0;
				uint index = 0;

				for (auto namedTag : compoundTag.tags)
				{
					if (namedTag.type == TAG_End) continue;

					if (namedTag.type == TAG_Int)
					{
						uint p1 = (uint)namedTag.payload[0];
						uint p2 = (uint)namedTag.payload[1] << 8;
						uint p3 = (uint)namedTag.payload[2] << 16;
						uint p4 = (uint)namedTag.payload[3] << 24;

						int pVal = p1 | p2 | p3 | p4;
						if (namedTag.name == "grid value") value = pVal;
						if (namedTag.name == "grid index") index = pVal;
					}
				}

				grid[index] = value;
			}
		}

		if (tag.name == "brick list")
		{
			for (auto compoundTag : tag.tags)
			{
				int brickBufferOffset = 0;

				for (auto namedTag : compoundTag.tags)
				{
					if (namedTag.type == TAG_Int)
					{
						uint payloadValue = (uint)namedTag.payload[0] | (uint)namedTag.payload[1] << 8 | (uint)namedTag.payload[2] << 16 | (uint)namedTag.payload[3] << 24;
						if (namedTag.name == "brick buffer offset") brickBufferOffset = payloadValue;
						if (namedTag.name == "brick zeroes") zeroes[brickBufferOffset]= payloadValue;
					}

					if (namedTag.type == TAG_Byte_Array)
					{
						if (namedTag.name == "brick value")
						{
							memcpy(brick + brickBufferOffset * BRICKSIZE, &namedTag.payload[0], BRICKSIZE * PAYLOADSIZE);
						}
					}
				}
				world.Mark(brickBufferOffset);
				world.NewBrick(); // Call New Brick to ensure trash buffer is initialised correctly
			}
		}
	}

	// Remove all states for a fresh Undo/Redo
	stateCurrent = stateHead;
	while (stateCurrent != stateTail)
	{
		stateTail = stateTail->prevState;
		DeleteState(stateTail->nextState);
		stateTail->nextState = NULL;
	}

	// Check for added lights in the scene
	vector<Light> vls;
	world.FindLightsInWord(vls);
	world.SetupLightBuffer(vls);

	// Optimize and forcve sync bricks to GPU
	world.OptimizeBricks();
	world.ForceSyncAllBricks();
}

void WorldEditor::SaveWorld()
{
	using namespace NBTHelper;

	World& world = *GetWorld();

	unsigned int* brick = world.GetBrick();
	uint* grid = world.GetGrid();
	uint* zeroes = world.GetZeroes();

	// NBT File format 
	//	Tag_List("grid list") : x entries of type Tag_Compound (Note: Elements in list are unnamed) 
	//	{
	//		TAG_Compound : 2 entries
	//		{
	//			TAG_Int("grid index") : index
	//			TAG_Int("grid value") : value
	//			TAG_End
	//		}	
	//		....
	//	}
	//  Tag_List("brick list") : x entries of type Tag_Compound 
	//	{
	//		TAG_Compound : 2 entries
	//		{
	//			TAG_Int(""brick buffer offset"") : offset
	//			TAG_Byte_Array("brick value") : [PAYLOADSIZE * BRICKSIZE bytes]
	// 			TAG_Int("brick zeroes") : zeroes
	//			TAG_End
	//		}	
	//		....
	//	}
	//   

	ofstream wf("worlddata.nbt", ios::out | ios::binary);
	if (!wf)
	{
		printf("Error opening file worlddata.nbt\n");
		return;
	}

	Tag listTag;
	listTag.type = TAG_List;
	listTag.name = "grid list";
	listTag.payload.push_back((byte)TAG_Compound);

	std::set<uint> oldBrickIdxs;
	std::unordered_map<uint, uint> oldToNew;
	uint newBrickIdx = 0;

	//	Loop through the world grid and see which indices/values need to be saved to file
	for (uint index = 0; index < GRIDHEIGHT * GRIDWIDTH * GRIDDEPTH; index++)
	{
		if (grid[index] == 0) continue;

		Tag compoundTag;
		compoundTag.type = TAG_Compound;
		compoundTag.name = "";

		Tag indexTag;
		indexTag.type = TAG_Int;
		indexTag.name = "grid index";
		union { byte bVal[4]; uint uVal; };
		uVal = index;
		indexTag.payload.assign(bVal, bVal + 4);

		Tag valueTag;
		valueTag.type = TAG_Int;
		valueTag.name = "grid value";

		// Rewire the brick offsets so they start at 0 and go to Num Unique Bricks
		if (!IsSolidGridCell(grid[index]))
		{
			int oldBrickOffset = grid[index] >> 1;
			oldBrickIdxs.insert(oldBrickOffset);
			auto newBrickVal = oldToNew.find(oldBrickOffset);
			if (newBrickVal == oldToNew.end())
			{
				oldToNew.insert(std::make_pair(oldBrickOffset, newBrickIdx));
				uVal = newBrickIdx++ << 1 | 1;
			}
			else
			{
				uVal = newBrickVal->second << 1 | 1;
			}

			valueTag.payload.assign(bVal, bVal + 4);
		}
		else
		{
			uVal = grid[index];
			valueTag.payload.assign(bVal, bVal + 4);

		}

		Tag endTag;
		endTag.type = TAG_End;

		compoundTag.tags.push_back(indexTag);
		compoundTag.tags.push_back(valueTag);
		compoundTag.tags.push_back(endTag);

		listTag.tags.push_back(compoundTag);
	}

	WriteTag(wf, listTag);

	listTag.tags.clear();
	listTag.name = "brick list";

	//	Loop through the unique set of bricks and save data to file
	for (auto oldBrickIdx : oldBrickIdxs)
	{
		Tag compoundTag;
		compoundTag.type = TAG_Compound;
		compoundTag.name = "";

		Tag offsetTag;
		offsetTag.type = TAG_Int;
		offsetTag.name = "brick buffer offset";
		union { byte bVal[4]; uint uVal; };
		uVal = oldToNew.find(oldBrickIdx)->second;
		offsetTag.payload.assign(bVal, bVal + 4);

		Tag valueTag;
		valueTag.type = TAG_Byte_Array;
		valueTag.name = "brick value";
		valueTag.payload.resize(PAYLOADSIZE * BRICKSIZE);
		memcpy(&valueTag.payload[0], brick + oldBrickIdx * BRICKSIZE, PAYLOADSIZE * BRICKSIZE);

		Tag zeroesTag;
		zeroesTag.type = TAG_Int;
		zeroesTag.name = "brick zeroes";
		uVal = zeroes[oldBrickIdx];
		zeroesTag.payload.assign(bVal, bVal + 4);

		NBTHelper::Tag endTag;
		endTag.type = TAG_End;

		compoundTag.tags.push_back(offsetTag);
		compoundTag.tags.push_back(valueTag);
		compoundTag.tags.push_back(zeroesTag);
		compoundTag.tags.push_back(endTag);

		listTag.tags.push_back(compoundTag);
	}
	NBTHelper::WriteTag(wf, listTag);

	wf.close();
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
#pragma endregion



#pragma region GUI
static void HelpMarker(const char* desc)
{
	ImGui::TextDisabled("(?)");
	if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
}

void WorldEditor::RenderGUI()
{
	if (!enabled) return;

	World& world = *GetWorld();
	RenderParams& params = world.GetRenderParams();
	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	static bool showCameraWindow = true;
	static bool showAllTiles = false;
	static bool drawGrid = true;

	// Saved indicies for proper ordering of the Tile Buttons
	static bool savedTileInit = true;
	static int savedTileIndices[8] = {};
	if (savedTileInit)
	{
		for (int n = 0; n < IM_ARRAYSIZE(savedTileIndices); n++)
		{
			savedTileIndices[n] = n;
		}
		savedTileInit = false;
	}

	// Saved indicies for proper ordering of the Big Tile Buttons
	static bool savedBigTileInit = true;
	static int savedBigTileIndices[8] = {};
	if (savedBigTileInit)
	{
		for (int n = 0; n < IM_ARRAYSIZE(savedBigTileIndices); n++)
		{
			savedBigTileIndices[n] = n;
		}
		savedBigTileInit = false;
	}

	// Saved indicies for proper ordering of the Sprite Buttons
	static bool savedSpriteInit = true;
	static int savedSpriteIndices[8] = {};
	if (savedSpriteInit)
	{
		for (int n = 0; n < IM_ARRAYSIZE(savedSpriteIndices); n++)
		{
			savedSpriteIndices[n] = n;
		}
		savedSpriteInit = false;
	}

	// Main Menu
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("New World")) { world.Clear(); world.Mark(0); }
			if (ImGui::MenuItem("Save World")) { SaveWorld(); }
			if (ImGui::MenuItem("Load World")) { LoadWorld(); } 

			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Edit"))
		{

			if (ImGui::MenuItem("Undo", "CTRL+Z", false, stateCurrent->prevState != NULL)) { Undo(); }
			if (ImGui::MenuItem("Redo", "CTRL+Y", false, stateCurrent->nextState != NULL)) { Redo(); }
			if (ImGui::MenuItem("Optimize Bricks", NULL, false)) { world.OptimizeBricks(); }
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("View"))
		{
			ImGui::MenuItem("Camera", NULL, &showCameraWindow);
			ImGui::MenuItem("Render Grid", NULL, &drawGrid);
			ImGui::MenuItem("Render Sprites", NULL, &drawGrid);
			ImGui::MenuItem("Show All Assets", NULL, &showAllTiles);
			ImGui::EndMenu();
		}

		if (ImGui::BeginMenu("Stats"))
		{
			ImGui::Text("Albedo render time: %.2f ms",  world.GetAlbedoTime() * 1000);
			ImGui::Text("Total render time: %.2f ms",  world.GetRenderTime() * 1000);
			ImGui::EndMenu();
		}

		ImGui::EndMainMenuBar();
	}

	// Edit Tool Window
	ImGui::Begin("Hot Bar");
	auto StyleHotBarTab = [](std::string tabName, std::vector<std::pair<int, GLuint>>& assets, int& assetIdx, int* savedIndex) -> bool
	{
		if (ImGui::BeginTabItem(tabName.c_str()))
		{
			int numButtons = min(8, (int)assets.size());
			for (int i = 0; i < numButtons; i++)
			{
				ImGui::PushID(i);
				int buttonIdx = savedIndex[i];
				ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(1.0, 0.0, 0.0, 1.0));
				if (buttonIdx == assetIdx)
					ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 5.0f);
				else
					ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 0.0f);


				if (ImGui::ImageButton((void*)(intptr_t)assets[buttonIdx].second, ImVec2(TILE_IMAGE_DIM / 2, TILE_IMAGE_DIM / 2), ImVec2(0.0f, 0.0f), ImVec2(1.0f, 1.0f)))
				{
					assetIdx = buttonIdx;
				}

				ImGui::PopStyleVar();
				ImGui::PopStyleColor();

				// Our buttons are both drag sources and drag targets here!
				if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None))
				{
					// Set payload to carry the index of our item (could be anything)
					ImGui::SetDragDropPayload("HotBar", &i, sizeof(int));
					ImGui::EndDragDropSource();
				}

				if (ImGui::BeginDragDropTarget())
				{
					if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("HotBar"))
					{
						IM_ASSERT(payload->DataSize == sizeof(int));
						int payload_n = *(const int*)payload->Data;
						int tmp = savedIndex[i];
						savedIndex[i] = savedIndex[payload_n];
						savedIndex[payload_n] = tmp;
					}

					if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(tabName.c_str()))
					{
						IM_ASSERT(payload->DataSize == sizeof(int));
						int payload_n = *(const int*)payload->Data;
						savedIndex[i] = payload_n;
						
					}
					ImGui::EndDragDropTarget();
				}

				ImGui::SameLine();
				ImGui::PopID();
			}
			// switch to the newly selected tab
			ImGui::EndTabItem();
			return true;
		}
		return false;
	};

	auto StyleVerticalTab = [](std::string tabName, std::vector<std::pair<int, GLuint>>& assets, int& assetIdx) -> bool
	{
		if (ImGui::BeginTabItem(tabName.c_str()))
		{
			for (int i = 0; i < (int)assets.size(); i++)
			{
				ImGui::PushID(i);

				// Center the Image Button for the vertical window
				ImGuiStyle& style = ImGui::GetStyle();
				float size = 128 + style.FramePadding.x * 2.0f;
				float avail = ImGui::GetContentRegionAvail().x;

				float off = (avail - size) * 0.5f;
				if (off > 0.0f)
					ImGui::SetCursorPosX(ImGui::GetCursorPosX() + off);

				ImGui::Image((void*)(intptr_t)assets[i].second, ImVec2(TILE_IMAGE_DIM / 2, TILE_IMAGE_DIM / 2), ImVec2(0.0f, 0.0f), ImVec2(1.0f, 1.0f));

				if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_SourceAllowNullID))
				{
					// Set payload to carry the index of our item (could be anything)
					ImGui::SetDragDropPayload(tabName.c_str(), &i, sizeof(int));
					ImGui::EndDragDropSource();
				}

				ImGui::PopID();
			}
			// switch to the newly selected tab
			ImGui::EndTabItem();
			return true;
		}
		return false;
	};

	if (ImGui::BeginTabBar("EditingTabBar", ImGuiTabBarFlags_AutoSelectNewTabs))
	{
		if (StyleHotBarTab("Tiles", loadedTiles, selectedTileIdx, savedTileIndices))
		{
			gesture.size = GestureSize::GESTURE_TILE;
		}
		if (StyleHotBarTab("Big Tiles", loadedBigTiles, selectedBigTileIdx, savedBigTileIndices))
		{
			gesture.size = GestureSize::GESTURE_BIG_TILE;
		}
		if (StyleHotBarTab("Sprites", loadedSprites, selectedSpriteIdx, savedSpriteIndices))
		{
			gesture.size = GestureSize::GESTURE_SPRITE;
		}

		if (ImGui::BeginTabItem("Brick/Voxel"))
		{

			static ImVec4 color = { 1.0f, 1.0f, 1.0f, 1.0f};

			static bool voxel = false;
			static bool refColor = false;
			static int displayMode = 0;
			static int pickerMode = 0;
			static int emissiveVal = 0;

			ImGuiColorEditFlags flags = ImGuiColorEditFlags_None;

			// Generate a default palette. The palette will persist and can be edited.
			static bool saved_palette_init = true;
			static ImVec4 saved_palette[32] = {};
			if (saved_palette_init)
			{
				for (int n = 0; n < IM_ARRAYSIZE(saved_palette); n++)
				{
					ImGui::ColorConvertHSVtoRGB(n / 31.0f, 0.8f, 0.8f,
						saved_palette[n].x, saved_palette[n].y, saved_palette[n].z);
					saved_palette[n].w = 1.0f; // Alpha
				}
				saved_palette_init = false;
			}

			ImGui::Columns(2);
			ImGui::SetColumnWidth(0, 150);
			static ImVec4 backup_color;
			ImGui::Text("Voxel Color:");
			bool picker = ImGui::ColorButton("##voxelcolor", color, flags, ImVec2(128, 128)); 
			ImGui::NextColumn();
			ImGui::NewLine();
			{
				ImGui::PushItemWidth(ImGui::GetColumnWidth(1) * 0.35);
				ImGui::SliderInt ("Emissive Strength", &emissiveVal, 0, 255, "%d", ImGuiSliderFlags_AlwaysClamp);
				ImGui::PopItemWidth();
			}
			ImGui::Checkbox("Edit Voxel", &voxel);
			ImGui::Columns();
			if (voxel) gesture.size = GestureSize::GESTURE_VOXEL;
			else gesture.size = GestureSize::GESTURE_BRICK;
			
			if (picker)
			{
				ImGui::OpenPopup("Popup Picker");
				backup_color = color;
			}

			if (ImGui::BeginPopup("Popup Picker"))
			{
				ImGui::Text("Brick Color Picker");
				ImGui::Separator();
				ImGui::ColorPicker4("##picker", (float*)&color, flags | ImGuiColorEditFlags_AlphaBar | ImGuiColorEditFlags_NoSmallPreview);
				ImGui::SameLine(); HelpMarker("Right-click on the individual color widget to show options.");
				ImGui::SameLine();

				ImGui::BeginGroup(); // Lock X position
				ImGui::Text("Current");
				ImGui::ColorButton("##current", color, ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_AlphaPreviewHalf, ImVec2(60, 40));
				ImGui::Text("Previous");
				if (ImGui::ColorButton("##previous", backup_color, ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_AlphaPreviewHalf, ImVec2(60, 40)))
					color = backup_color;
				ImGui::Separator();
				ImGui::Text("Palette"); ImGui::SameLine(); HelpMarker("Drag and drop colors from the picker to update the palette.");
				for (int n = 0; n < IM_ARRAYSIZE(saved_palette); n++)
				{
					ImGui::PushID(n);
					if ((n % 8) != 0)
						ImGui::SameLine(0.0f, ImGui::GetStyle().ItemSpacing.y);

					ImGuiColorEditFlags palette_button_flags = ImGuiColorEditFlags_NoAlpha | ImGuiColorEditFlags_NoPicker | ImGuiColorEditFlags_NoTooltip;
					if (ImGui::ColorButton("##palette", saved_palette[n], palette_button_flags, ImVec2(20, 20)))
						color = ImVec4(saved_palette[n].x, saved_palette[n].y, saved_palette[n].z, saved_palette[n].w);

					// Allow user to drop colors into each palette entry. Note that ColorButton() is already a
					// drag source by default, unless specifying the ImGuiColorEditFlags_NoDragDrop flag.
					if (ImGui::BeginDragDropTarget())
					{
						if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(IMGUI_PAYLOAD_TYPE_COLOR_3F))
							memcpy((float*)&saved_palette[n], payload->Data, sizeof(float) * 3);
						if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(IMGUI_PAYLOAD_TYPE_COLOR_4F))
							memcpy((float*)&saved_palette[n], payload->Data, sizeof(float) * 4);
						ImGui::EndDragDropTarget();
					}

					ImGui::PopID();
				}
				ImGui::EndGroup();
				ImGui::EndPopup();
			}

			const uint red = min(15u, (uint)(color.x * 15.0f));
			const uint green = min(15u, (uint)(color.y * 15.0f));
			const uint blue = min(15u, (uint)(color.z * 15.0f));
			const uint alpha = min(15u, (uint)(color.w * 15.0f));
			const uint emissiveness = min(255u, (uint)(emissiveVal));

			color.x = (red * (1.0f / 15.0f));
			color.y = (green * (1.0f / 15.0f));
			color.z = (blue * (1.0f / 15.0f));
			color.w = (alpha * (1.0f / 15.0f));

			const uint p = (red << 8) + (green << 4) + blue;
			voxelValue = (p == 0) ? 1 : p;

			voxelValue |= alpha << 12;
			voxelValue |= emissiveness << 16;
			
			ImGui::EndTabItem();
		}

		ImGui::EndTabBar();
	}
	ImGui::End();

	// Camera Window
	if (showCameraWindow)
	{
		ImGui::Begin("Camera");
		float cameraPos[3];
		cameraPos[0] = world.GetCameraPos().x;	
		cameraPos[1] = world.GetCameraPos().y;
		cameraPos[2] = world.GetCameraPos().z;

		ImGui::InputFloat3("Position", cameraPos, "%.1f", ImGuiInputTextFlags_CharsDecimal | ImGuiInputTextFlags_EnterReturnsTrue);
		world.SetCameraPos(make_float3(cameraPos[0], cameraPos[1], cameraPos[2]));


		float cameraViewDir[3];
		cameraViewDir[0] = world.GetCameraViewDir().x;
		cameraViewDir[1] = world.GetCameraViewDir().y;
		cameraViewDir[2] = world.GetCameraViewDir().z;
		ImGui::InputFloat3("Direction", cameraViewDir, "%.3f", ImGuiInputTextFlags_CharsDecimal | ImGuiInputTextFlags_EnterReturnsTrue);
		world.SetCameraViewDir(make_float3(cameraViewDir[0], cameraViewDir[1], cameraViewDir[2]));
		ImGui::End();
	}

	if (showAllTiles)
	{
		ImGui::Begin("All Assets"); ImGui::SameLine(); HelpMarker("Drag and drop assets from here into the Hot Bar for use");
		if (ImGui::BeginTabBar("AllAssets", ImGuiTabBarFlags_AutoSelectNewTabs))
		{
			if (StyleVerticalTab("Tiles", loadedTiles, selectedTileIdx)) { gesture.size == GestureSize::GESTURE_TILE; }
			if (StyleVerticalTab("Big Tiles", loadedBigTiles, selectedBigTileIdx)) { gesture.size == GestureSize::GESTURE_BIG_TILE; }
			if (StyleVerticalTab("Sprites", loadedSprites, selectedSpriteIdx)) { gesture.size == GestureSize::GESTURE_SPRITE;  }
		}
		ImGui::EndTabBar();
		ImGui::End();
	}

	params.drawGrid = drawGrid;

	// Render dear imgui into screen
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
#pragma endregion