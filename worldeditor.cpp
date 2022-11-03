#include "precomp.h"
#include "worldeditor.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "lib/stb_image_write.h"

#define MAX_ALLOCATION 512000000 // max bytes allowed for undo/redo states

#define BIG_TILE_IMAGE_SCALE 8
#define TILE_IMAGE_SCALE 16

#define TILE_IMAGE_DIM BRICKDIM * TILE_IMAGE_SCALE

#pragma region NBTHelper
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
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); // This is required on WebGL for non power-of-two textures
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE); // Same

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

void WorldEditor::LoadTiles()
{
	const std::filesystem::path tiles{ "assets/Tiles" };
	const std::filesystem::path bigTiles{ "assets/BigTiles" };

	TileManager* tileManager = TileManager::GetTileManager();

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
				int imageSize = TILE_IMAGE_DIM * TILE_IMAGE_DIM;
				uint32_t* buffer = new uint32_t[imageSize];
				memset(buffer, 0, imageSize * sizeof(uint32_t));

				PAYLOAD* voxels = tileManager->tile[tileIdx]->voxels;

				// Find the first non transparent voxel for each x,y position
				for (int y = 0; y < BRICKDIM; y++)
				{
					for (int x = 0; x < BRICKDIM; x++)
					{
						for (int z = 0; z < BRICKDIM; z++)
						{
							PAYLOAD v = voxels[x + y * BRICKDIM + z * BRICKDIM * BRICKDIM];
							if (v)
							{
								float3 val = ToFloatRGB(v);

								// Fill in a square of pixels (TILE_IMAGE_SCALE * TILE_IMAGE_SCALE) for a given voxel color 
								for (int _x = x * TILE_IMAGE_SCALE; _x < (x + 1) * TILE_IMAGE_SCALE; _x++)
									for (int _y = y * TILE_IMAGE_SCALE; _y < (y + 1) * TILE_IMAGE_SCALE; _y++)
									{
										uint32_t& dst = buffer[((TILE_IMAGE_DIM - 1) - (_y)) * TILE_IMAGE_DIM + (TILE_IMAGE_DIM - 1) - (_x)];
										dst = ((std::min<uint32_t>(static_cast<uint32_t>(val.x * 255), 255) << 0) |
											(std::min<uint32_t>(static_cast<uint32_t>(val.y * 255), 255) << 8) |
											(std::min<uint32_t>(static_cast<uint32_t>(val.z * 255), 255) << 16) | 255 << 24);
									}

								break;
							}
						}

					}
				}

				stbi_write_png(p.string().c_str(), TILE_IMAGE_DIM, TILE_IMAGE_DIM, 4, buffer, TILE_IMAGE_DIM * sizeof(uint32_t));
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
				int imageSize = TILE_IMAGE_DIM * TILE_IMAGE_DIM;
				uint32_t* buffer = new uint32_t[imageSize];
				memset(buffer, 0, imageSize * sizeof(uint32_t));

				Tile* tiles = tileManager->bigTile[bigTileIdx]->tile;

				for (int subTile = 0; subTile < 8; subTile++)
				{
					int sx = subTile & 1, sy = (subTile >> 1) & 1, sz = (subTile >> 2) & 1;

					for (int z = BRICKDIM - 1; z > 0; z--) for (int y = 0; y < BRICKDIM; y++) for (int x = 0; x < BRICKDIM; x++)
					{
						PAYLOAD v = tiles[subTile].voxels[x + y * BRICKDIM + z * BRICKDIM * BRICKDIM];
						if (v)
						{
							float3 val = ToFloatRGB(v);

							int xPos = x * BIG_TILE_IMAGE_SCALE + sx * BIG_TILE_IMAGE_SCALE * BRICKDIM;
							int yPos = y * BIG_TILE_IMAGE_SCALE + sy * BIG_TILE_IMAGE_SCALE * BRICKDIM;
							// Fill in a square of pixels (BIG_TILE_IMAGE_SCALE * BIG_TILE_IMAGE_SCALE) for a given voxel color 
							for (int _x = xPos; _x < xPos + BIG_TILE_IMAGE_SCALE; _x++)
								for (int _y = yPos; _y < yPos + BIG_TILE_IMAGE_SCALE; _y++)
								{
									uint32_t& dst = buffer[((TILE_IMAGE_DIM - 1) - _y) * TILE_IMAGE_DIM + (TILE_IMAGE_DIM - 1) - _x];
									dst = ((std::min<uint32_t>(static_cast<uint32_t>(val.x * 255), 255) << 0) |
										(std::min<uint32_t>(static_cast<uint32_t>(val.y * 255), 255) << 8) |
										(std::min<uint32_t>(static_cast<uint32_t>(val.z * 255), 255) << 16) | 255 << 24);
								}
						}
					}
				}
				stbi_write_png(p.string().c_str(), TILE_IMAGE_DIM, TILE_IMAGE_DIM, 4, buffer, TILE_IMAGE_DIM * sizeof(uint32_t));
			}

			int my_image_width = 0;
			int my_image_height = 0;
			GLuint my_image_texture = 0;

			bool ret = LoadTextureFromFile(p.string().c_str(), &my_image_texture, &my_image_width, &my_image_height);
			IM_ASSERT(ret);
			loadedBigTiles.push_back(std::make_pair(bigTileIdx, my_image_texture));
		}
	}
}

WorldEditor::WorldEditor()
{
	tempBricks = (PAYLOAD*)_aligned_malloc(CHUNKCOUNT * CHUNKSIZE, 64);
	tempGrid = (uint*)_aligned_malloc(GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * 4, 64);
	//tempBrickInfo = (BrickInfo*)_aligned_malloc(BRICKCOUNT * sizeof(BrickInfo), 64);
	stateHead = (State*)calloc(1, sizeof(State));
	stateCurrent = stateTail = stateHead;

	memset(tempBricks, 0, CHUNKCOUNT * CHUNKSIZE);
	memset(tempGrid, 0, GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));
	//memset(tempBrickInfo, BRICKSIZE, BRICKCOUNT * sizeof(BrickInfo));

	LoadTiles();
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

	if (gesture.mode & GestureMode::GESTURE_MULTI)
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
		//memcpy(tempBrickInfo, world.GetBrickInfo(), GRIDWIDTH * GRIDHEIGHT * GRIDDEPTH * sizeof(uint));

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
void WorldEditor::UpdateEditedBricks(int x, int y, int z)
{
	if (gesture.size == GestureSize::GESTURE_BIG_TILE)
	{
		for (int _x = 0; _x < 2; _x++)
			for (int _y = 0; _y < 2; _y++)
				for (int _z = 0; _z < 2; _z++)
					editedBricks.insert(make_int3(x / BRICKDIM + _x, y / BRICKDIM + _y, z / BRICKDIM + _z));
	}
	else
	{
		editedBricks.insert((make_int3)(x / BRICKDIM, y / BRICKDIM, z / BRICKDIM));
	}

	return;
}

// Allow the adding/removing of multiple bricks in the seclected box
void WorldEditor::MultiAddRemove()
{
	aabb oldBox = selected.box;
	selected.box = selected.anchor;
	UpdateSelectedBox();

	if (oldBox == selected.box) return;

	World& world = *GetWorld();
	unsigned int* brick = world.GetBrick();
	uint* grid = world.GetGrid();

	aabb newBox = selected.box;
	// Compute the overlap betweeen the old and new aabbs

	int boxScale = 1;
	switch (gesture.size)
	{
	case GestureSize::GESTURE_VOXEL:
		boxScale = 1;
		break;
	case GestureSize::GESTURE_TILE:
	case GestureSize::GESTURE_BRICK:
		boxScale = BRICKDIM;
		break;
	case GestureSize::GESTURE_BIG_TILE:
		boxScale = BRICKDIM * 2;
		break;
	default:
		break;
	}


	oldBox.bmin3 = make_float3((int)oldBox.bmin3.x / boxScale, (int)oldBox.bmin3.y / boxScale, (int)oldBox.bmin3.z / boxScale);
	oldBox.bmax3 = make_float3((int)oldBox.bmax3.x / boxScale, (int)oldBox.bmax3.y / boxScale, (int)oldBox.bmax3.z / boxScale);
	newBox.bmin3 = make_float3((int)newBox.bmin3.x / boxScale, (int)newBox.bmin3.y / boxScale, (int)newBox.bmin3.z / boxScale);
	newBox.bmax3 = make_float3((int)newBox.bmax3.x / boxScale, (int)newBox.bmax3.y / boxScale, (int)newBox.bmax3.z / boxScale);
	aabb intersection = oldBox.Intersection(newBox);

	auto RestoreTempBrickVal = [&](int x, int y, int z) {
		const uint cellIdx = x + z * GRIDWIDTH + y * GRIDWIDTH * GRIDDEPTH;
		const uint tempGridVal = tempGrid[cellIdx];
		const uint curGridVal = grid[cellIdx];

		// If we're reintroducing what was an empty/solid brick, just remove the current one
		if ((tempGridVal & 1) == 0)
		{
			world.RemoveBrick(x, y, z);
			grid[cellIdx] = tempGridVal;
			return;
		}

		// Get current and new brick offsets (g1)
		uint tempBrickOffset = (tempGridVal >> 1);
		uint curBrickOffset = (curGridVal >> 1);

		uint brickIdx = tempBrickOffset;
		uint gridValue = tempGridVal;

		// If the current brick is empty/solid we need to create a NewBrick
		if ((curGridVal & 1) == 0) brickIdx = world.NewBrick(), gridValue = (brickIdx << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(brick + brickIdx * BRICKSIZE, tempBricks + tempBrickOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
		grid[cellIdx] = gridValue;
		//world.GetBrickInfo()[brickIdx].zeroes = tempBrickInfo[tempBrickOffset].zeroes;
		world.Mark(brickIdx);
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
					const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
					// obtain current brick identifier from top-level grid
					uint g = grid[cellIdx], g1 = g >> 1;

					// calculate the position of the voxel inside the brick
					const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
					const uint voxelIdx = g1 * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM;
					const uint oldVoxelValue = tempBricks[voxelIdx];
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
					UpdateEditedBricks(x * boxScale, y * boxScale, z * boxScale);
					continue;
				}

				if (gesture.mode & GestureMode::GESTURE_REMOVE)
				{
					Remove(x * boxScale, y * boxScale, z * boxScale);
				}
				else
				{
					Add(x * boxScale, y * boxScale, z * boxScale);
				}
			}
}

void WorldEditor::Add(int vx, int vy, int vz)
{
	World& world = *GetWorld();

	// Convert voxel xyz to brick coords
	int bx = vx / BRICKDIM;
	int by = vy / BRICKDIM;
	int bz = vz / BRICKDIM;

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
		}
	}
}

void WorldEditor::Remove(int vx, int vy, int vz)
{
	World& world = *GetWorld();

	int bx = vx / BRICKDIM;
	int by = vy / BRICKDIM;
	int bz = vz / BRICKDIM;

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
		}
	}

}

// Update which brick is currently selected by the mouse cursor 
void WorldEditor::UpdateSelectedBox()
{
	World& world = *GetWorld();
	RenderParams& params = world.GetRenderParams();

	int boxScale = 1;
	switch (gesture.size)
	{
	case GestureSize::GESTURE_VOXEL:
		boxScale = 1;
		break;
	case GestureSize::GESTURE_TILE:
	case GestureSize::GESTURE_BRICK:
		boxScale = BRICKDIM;
		break;
	case GestureSize::GESTURE_BIG_TILE:
		boxScale = BRICKDIM * 2;
		break;
	default:
		break;
	}

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


		float3 hitPoint = ray.O + ray.D * distance;
		// Get position inside of the grid to determine brick location
		float3 gridPos = hitPoint + 0.5 * normal;

		if (gridPos.x < 0 || gridPos.y < 0 || gridPos.z < 0 || gridPos.x > MAPWIDTH || gridPos.y > MAPHEIGHT || gridPos.z > MAPDEPTH)
			return;

		float3 boxPos;
		boxPos = make_float3(
			(int)(gridPos.x / boxScale) * boxScale,
			(int)(gridPos.y / boxScale) * boxScale,
			(int)(gridPos.z / boxScale) * boxScale
		);
		selected.box.Grow(boxPos);
		selected.box.Grow(boxPos + make_float3(boxScale - 1, boxScale - 1, boxScale - 1));
	}
	else
	{
		float t = intersection.GetDistance();
		float3 normal = intersection.GetNormal();
		float3 hitPoint = ray.O + ray.D * t;

		// Get position inside of the voxel to determine brick location
		float3 voxelPos = hitPoint - 0.1 * normal;
		float3 boxPos;
		boxPos = make_float3(
			(int)(voxelPos.x / boxScale) * boxScale,
			(int)(voxelPos.y / boxScale) * boxScale,
			(int)(voxelPos.z / boxScale) * boxScale
		);

		if (gesture.mode & GestureMode::GESTURE_REMOVE)
		{
			selected.box.Grow(boxPos);
			selected.box.Grow(boxPos + make_float3(boxScale - 1, boxScale - 1, boxScale - 1));

		}
		else
		{
			float3 newBoxPos = boxPos + normal * boxScale;

			// If the normal puts us outside of the grid, just use the original intersection
			if (newBoxPos.x < 0 || newBoxPos.y < 0 || newBoxPos.z < 0 || newBoxPos.x > MAPWIDTH || newBoxPos.y > MAPHEIGHT || newBoxPos.z > MAPDEPTH)
			{
				selected.box.Grow(boxPos);
				selected.box.Grow(boxPos + make_float3(boxScale - 1, boxScale - 1, boxScale - 1));
			}
			else
			{
				selected.box.Grow(newBoxPos);
				selected.box.Grow(newBoxPos + make_float3(boxScale - 1, boxScale - 1, boxScale - 1));
			}
		}
	}

	// Update rendering params to trace the selected box outline
	params.selectedMin = selected.box.bmin3;
	params.selectedMax = selected.box.bmax3;
	if (gesture.size == GestureSize::GESTURE_BIG_TILE) params.wireBoxWidth = 0.5f;
	if (gesture.size == GestureSize::GESTURE_TILE || GESTURE_BRICK) params.wireBoxWidth = 0.3f;
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
	//BrickInfo* brickInfo = world.GetBrickInfo();

	int numBricks = stateCurrent->numBricks;
	for (int idx = 0; idx < numBricks; idx++)
	{
		int3 b = stateCurrent->editedBricks[idx];
		const uint cellIdx = b.x + b.z * GRIDWIDTH + b.y * GRIDWIDTH * GRIDDEPTH;

		const uint oldGridVal = stateCurrent->oldGridVals[idx];
		const uint newGridVal = stateCurrent->newGridVals[idx];

		// If we're restoring what was an empty/solid brick, just remove the current one
		if ((oldGridVal & 1) == 0)
		{
			world.RemoveBrick(b.x, b.y, b.z);
			grid[cellIdx] = oldGridVal;
			continue;
		}

		// Get old and new brick offsets (g1)
		uint oldBrickOffset = (oldGridVal >> 1);
		uint newBrickOffset = (newGridVal >> 1);

		uint brickIdx = oldBrickOffset;
		uint gridValue = oldGridVal;

		// If the current brick is empty/solid we need to create a NewBrick
		if ((newGridVal & 1) == 0) brickIdx = world.NewBrick(), gridValue = (brickIdx << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(bricks + (brickIdx * BRICKSIZE), stateCurrent->oldBricks + (idx * BRICKSIZE), BRICKSIZE * PAYLOADSIZE);
		grid[cellIdx] = gridValue;
		//brickInfo[brickIdx].zeroes = stateCurrent->oldBrickZeroes[idx];
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
	//BrickInfo* brickInfo = world.GetBrickInfo();

	int numBricks = stateCurrent->numBricks;
	for (int idx = 0; idx < numBricks; idx++)
	{
		int3 b = stateCurrent->editedBricks[idx];
		const uint cellIdx = b.x + b.z * GRIDWIDTH + b.y * GRIDWIDTH * GRIDDEPTH;

		const uint oldGridVal = stateCurrent->oldGridVals[idx];
		const uint newGridVal = stateCurrent->newGridVals[idx];

		// If we're restoring what was an empty/solid brick, just remove the current one
		if ((newGridVal & 1) == 0)
		{
			world.RemoveBrick(b.x, b.y, b.z);
			grid[cellIdx] = newGridVal;
			continue;
		}

		// Get old and new brick offsets (g1)
		uint oldBrickOffset = (oldGridVal >> 1);
		uint newBrickOffset = (newGridVal >> 1);

		uint brickIdx = newBrickOffset;
		uint gridValue = newGridVal;

		// If the current brick is empty/solid we need to create a NewBrick
		if ((oldGridVal & 1) == 0) brickIdx = world.NewBrick(), gridValue = (brickIdx << 1) | 1;

		// Copy the brick from the saved temporary brick buffer to our current brick buffer 
		memcpy(bricks + (brickIdx * BRICKSIZE), stateCurrent->newBricks + (idx * BRICKSIZE), BRICKSIZE * PAYLOADSIZE);
		grid[cellIdx] = gridValue;
		//brickInfo[brickIdx].zeroes = stateCurrent->oldBrickZeroes[idx];
		world.Mark(brickIdx);
	}

	//vector<Light> ls;
	//world.SetupLights(ls);
}

// Called after a gesture has been completed so we can update the history
void WorldEditor::SaveState()
{
	World& world = *GetWorld();
	uint* grid = world.GetGrid();
	PAYLOAD* bricks = world.GetBrick();
	//BrickInfo* brickInfo = world.GetBrickInfo();
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
	stateCurrent->newGridVals = (uint*)_aligned_malloc(numBricks * 4, 64);
	stateCurrent->oldGridVals = (uint*)_aligned_malloc(numBricks * 4, 64);
	//stateCurrent->newBrickZeroes = (uint*)_aligned_malloc(numBricks * sizeof(uint), 64);
	//stateCurrent->oldBrickZeroes = (uint*)_aligned_malloc(numBricks * sizeof(uint), 64);
	stateCurrent->editedBricks = (int3*)_aligned_malloc(numBricks * sizeof(int3), 64);

	int idx = 0;
	for (auto itr = editedBricks.begin(); itr != editedBricks.end(); itr++)
	{
		int3 brick = *itr;
		const uint cellIdx = brick.x + brick.z * GRIDWIDTH + brick.y * GRIDWIDTH * GRIDDEPTH;

		const uint tempGridValue = tempGrid[cellIdx];
		const uint gridValue = grid[cellIdx];

		uint tempBrickOffset = (tempGridValue >> 1);
		uint brickOffset = (gridValue >> 1);

		// Old brick is empty/or a solid brick
		if ((tempGridValue & 1) == 0)
		{
			memset(stateCurrent->oldBricks + (idx * BRICKSIZE), 0, BRICKSIZE * PAYLOADSIZE);
		}
		else
		{
			memcpy(stateCurrent->oldBricks + (idx * BRICKSIZE), tempBricks + tempBrickOffset * BRICKSIZE, BRICKSIZE * PAYLOADSIZE);
		}
		
		// New brick is empty/or a solid brick
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

		//stateCurrent->oldBrickZeroes[idx] = tempBrickInfo[tempBrickOffset].zeroes;
		//stateCurrent->newBrickZeroes[idx] = brickInfo[brickOffset].zeroes;

		stateCurrent->editedBricks[idx] = brick;
		idx++;
	}


	// Add the amount of allocated memory 
	allocatedUndo += 2 * (numBricks * BRICKSIZE * PAYLOADSIZE);
	allocatedUndo += 2 * (numBricks * 4);
	allocatedUndo += 2 * (numBricks * sizeof(uint));
	allocatedUndo += (numBricks * sizeof(int3));


	// Check to see if the state has actually changed....if not don't add it to the history
	if (memcmp(stateCurrent->oldGridVals, stateCurrent->newGridVals, numBricks * 4) == 0)
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
	//BrickInfo* brickInfo = world.GetBrickInfo();

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
				int brickIndex = 0;

				for (auto namedTag : compoundTag.tags)
				{
					if (namedTag.type == TAG_Int)
					{
						uint payloadValue = (uint)namedTag.payload[0] | (uint)namedTag.payload[1] << 8 | (uint)namedTag.payload[2] << 16 | (uint)namedTag.payload[3] << 24;
						if (namedTag.name == "brick index") brickIndex = payloadValue;
						//if (namedTag.name == "brick zeroes") brickInfo[brickIndex].zeroes = payloadValue;
					}

					if (namedTag.type == TAG_Byte_Array)
					{
						if (namedTag.name == "brick value")
						{
							memcpy(brick + brickIndex * BRICKSIZE, &namedTag.payload[0], BRICKSIZE * PAYLOADSIZE);
						}
					}
				}
				world.Mark(brickIndex); // Mark to send to GPU
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
}

void WorldEditor::SaveWorld()
{
	using namespace NBTHelper;

	World& world = *GetWorld();

	unsigned int* brick = world.GetBrick();
	uint* grid = world.GetGrid();
	//BrickInfo* brickInfo = world.GetBrickInfo();

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
	//			TAG_Int("brick index") : index
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
	for (uint idx = 0; idx < GRIDHEIGHT * GRIDWIDTH * GRIDDEPTH; idx++)
	{
		if (grid[idx] == 0) continue;

		Tag compoundTag;
		compoundTag.type = TAG_Compound;
		compoundTag.name = "";

		Tag indexTag;
		indexTag.type = TAG_Int;
		indexTag.name = "grid index";
		union { byte bVal[4]; uint uVal; };
		uVal = idx;
		indexTag.payload.assign(bVal, bVal + 4);

		Tag valueTag;
		valueTag.type = TAG_Int;
		valueTag.name = "grid value";

		// Rewire the brick indices so they start at 0 and go to Num Unique Bricks
		if (grid[idx] & 1)
		{
			int oldBrickIdx = grid[idx] >> 1;
			oldBrickIdxs.insert(oldBrickIdx);
			auto newBrickVal = oldToNew.find(oldBrickIdx);
			if (newBrickVal == oldToNew.end())
			{
				oldToNew.insert(std::make_pair(oldBrickIdx, newBrickIdx));
				uVal = newBrickIdx++ << 1 | 1;
			}
			else
			{
				uVal = newBrickVal->second << 1 | 1;
			}

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

		Tag indexTag;
		indexTag.type = TAG_Int;
		indexTag.name = "brick index";
		union { byte bVal[4]; uint uVal; };
		uVal = oldToNew.find(oldBrickIdx)->second;
		indexTag.payload.assign(bVal, bVal + 4);

		Tag valueTag;
		valueTag.type = TAG_Byte_Array;
		valueTag.name = "brick value";
		valueTag.payload.resize(PAYLOADSIZE * BRICKSIZE);
		memcpy(&valueTag.payload[0], brick + oldBrickIdx * BRICKSIZE, PAYLOADSIZE * BRICKSIZE);

		//Tag zeroesTag;
		//zeroesTag.type = TAG_Int;
		//zeroesTag.name = "brick zeroes";
		//uVal = brickInfo[oldBrickIdx].zeroes;
		//zeroesTag.payload.assign(bVal, bVal + 4);

		NBTHelper::Tag endTag;
		endTag.type = TAG_End;

		compoundTag.tags.push_back(indexTag);
		compoundTag.tags.push_back(valueTag);
		//compoundTag.tags.push_back(zeroesTag);
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

	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Save World")) { SaveWorld(); }
			if (ImGui::MenuItem("Load World")) { LoadWorld(); } 

			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Edit"))
		{

			if (ImGui::MenuItem("Undo", "CTRL+Z", false, stateCurrent->prevState != NULL)) { Undo(); }
			if (ImGui::MenuItem("Redo", "CTRL+Y", false, stateCurrent->nextState != NULL)) { Redo(); }
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}

	// render your GUI
	ImGui::Begin("Edit Tool");

	auto StyleTileTab = [](std::string tabName, std::vector<std::pair<int, GLuint>>& tiles, int& tileIdx) -> bool
	{
		if (ImGui::BeginTabItem(tabName.c_str()))
		{
			int numButtons = min(10, (int)tiles.size());
			for (int buttonIdx = 0; buttonIdx < numButtons; buttonIdx++)
			{
				ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(1.0, 0.0, 0.0, 1.0));
				if (buttonIdx == tileIdx)
					ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 5.0f);
				else
					ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 0.0f);


				if (ImGui::ImageButton((void*)(intptr_t)tiles[buttonIdx].second, ImVec2(TILE_IMAGE_DIM, TILE_IMAGE_DIM), ImVec2(0.0f, 0.0f), ImVec2(1.0f, 1.0f)))
				{
					tileIdx = buttonIdx;
				}

				ImGui::PopStyleVar();
				ImGui::PopStyleColor();
				ImGui::SameLine();
			}
			// switch to the newly selected tab
			ImGui::EndTabItem();
			return true;
		}
		return false;
	};

	ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_AutoSelectNewTabs;
	if (ImGui::BeginTabBar("EditingTabBar", tab_bar_flags))
	{
		if (StyleTileTab("Tile", loadedTiles, selectedTileIdx))
		{
			gesture.size = GestureSize::GESTURE_TILE;
		}
		if (StyleTileTab("Big Tile", loadedBigTiles, selectedBigTileIdx))
		{
			gesture.size = GestureSize::GESTURE_BIG_TILE;
		}

		if (ImGui::BeginTabItem("Brick/Voxel"))
		{

			static ImVec4 color = { 1.0f, 1.0f, 1.0f, 1.0f};

			static bool emmisive = false;
			static bool voxel = false;
			static bool ref_color = false;
			static int display_mode = 0;
			static int picker_mode = 0;


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

			static ImVec4 backup_color;
			ImGui::Text("Voxel Color:");
			bool picker = ImGui::ColorButton("##voxelcolor", color, flags, ImVec2(128, 128));
			ImGui::SameLine();

			ImGui::Checkbox("Emmisive", &emmisive);
			ImGui::Checkbox("Edit Voxel", &voxel);

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
						color = ImVec4(saved_palette[n].x, saved_palette[n].y, saved_palette[n].z, color.w); // Preserve alpha!

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

			color.x = (red * (1.0f / 15.0f));
			color.y = (green * (1.0f / 15.0f));
			color.z = (blue * (1.0f / 15.0f));
			color.w = (alpha * (1.0f / 15.0f));

			const uint p = (red << 8) + (green << 4) + blue;
			voxelValue = (p == 0) ? 1 : p;

			if (emmisive) voxelValue |= 15 << 12;
			ImGui::EndTabItem();
		}

		ImGui::EndTabBar();
	}

	ImGui::End();
	// Render dear imgui into screen
	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
#pragma endregion