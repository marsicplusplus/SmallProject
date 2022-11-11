#pragma once

namespace Tmpl8
{
	namespace NBTHelper
	{
		enum TagType {
			TAG_End = 0,
			TAG_Byte = 1,
			TAG_Short = 2,
			TAG_Int = 3,
			TAG_Long = 4,
			TAG_Float = 5,
			TAG_Double = 6,
			TAG_Byte_Array = 7,
			TAG_String = 8,
			TAG_List = 9,
			TAG_Compound = 10
		};

		struct Tag {
			int type = TAG_End;
			std::string name;
			std::vector<Tag> tags;
			std::vector<byte> payload;
		};

		void WriteTag(std::ofstream& wf, NBTHelper::Tag& tag);
		void WriteTagList(std::ofstream& wf, NBTHelper::Tag& tag);
		void WriteTagCompound(std::ofstream& wf, Tag& tag);
		void WriteTagType(std::ofstream& wf, int tagType);
		void WriteTagInt(std::ofstream& wf, NBTHelper::Tag& tag);
		void WriteTagName(std::ofstream& wf, std::string value);
		void WriteTagByteArray(std::ofstream& wf, Tag& tag);
		void ReadTagName(std::ifstream& rf, string& tagName);
		void ReadTagType(std::ifstream& rf, int& tagType);
		void ReadTag(std::ifstream& wf, Tag& tag, bool readTagName = true);
		void ReadTagInt(std::ifstream& rf, Tag& tag);
		void ReadTagList(std::ifstream& rf, Tag& tag);
		void ReadTagByteArray(std::ifstream& rf, Tag& tag);
		void ReadTagCompound(std::ifstream& rf, Tag& tag);

	}

	class WorldEditor
	{

		typedef struct State;
		struct State {
			PAYLOAD* newBricks;
			PAYLOAD* oldBricks;

			uint* oldCellValues;
			uint* newCellValues;

			uint* newBrickZeroes;
			uint* oldBrickZeroes;

			State* prevState = NULL;
			State* nextState = NULL;

			int3* editedBricks;
			uint numBricks;
		};

		struct Selected
		{
			aabb box;
			aabb anchor;
		};

		// Taking inspiration from Goxel's gesture class
		enum GestureButton {      
			GESTURE_NO_BUTTONS = 0,
			GESTURE_LMB = 1,
			GESTURE_RMB = 2
		};

		enum GestureKey {
			GESTURE_NO_KEYS = 0,
			GESTURE_CTRL = 1 << 0,
			GESTURE_SHIFT = 1 << 1
		};

		enum GestureState {
			GESTURE_POSSIBLE = 0,
			GESTURE_ACTIVE
		};

		enum GestureMode {
			GESTURE_ADD = 0,
			GESTURE_REMOVE = 1 << 0,
			GESTURE_MULTI = 1 << 1
		};

		enum GestureSize {
			GESTURE_VOXEL = 0,
			GESTURE_BRICK = 1,
			GESTURE_TILE = 2,
			GESTURE_BIG_TILE = 3,
			GESTURE_SPRITE = 4
		};

		struct Gesture
		{
			uint buttons = GestureButton::GESTURE_NO_BUTTONS;
			uint keys = GestureKey::GESTURE_NO_KEYS;
			uint state = GestureState::GESTURE_POSSIBLE;
			uint mode = GestureMode::GESTURE_ADD;
			uint size = GestureSize::GESTURE_TILE;
		};

	public:
		WorldEditor();
		void MouseMove(int x, int y);
		void MouseDown(int button);
		void MouseUp(int button);
		void KeyUp(int key);
		void KeyDown(int key);
		bool IsEnabled() { return enabled; }
		void Enable() { UpdateSelectedBox(); enabled = true; }
		void Disable() { ResetEditor(); enabled = false; }
		void ResetEditor();
		void RenderGUI();

	private:
		void UpdateSelectedBox();
		void Add(uint x, uint y, uint z);
		void Remove(uint x, uint y, uint z);
		void MultiAddRemove();
		void UpdateGestureMode();
		void CheckMemoryAllowance();
		void Redo();
		void Undo();
		void SaveState();
		void DeleteState(State* state);
		bool CreateNewState();
		void UpdateEditedBricks(uint bx, uint by, uint bz);
		void SaveWorld();
		void LoadWorld();
		int3 GetBoxScale();
		void AddBackLights(uint bx, uint by, uint bz);
		void LoadAssets();

		// Input and Gesture 
		int2 mousePos;
		Gesture gesture;
		uint selectedButtons = GestureButton::GESTURE_NO_BUTTONS;
		uint selectedKeys = GestureKey::GESTURE_NO_KEYS;
		Selected selected;

		// Temporary buffers to hold previous state
		PAYLOAD* tempBricks = 0;
		uint* tempGrid = 0;
		uint* tempZeroes = 0;

		std::vector<std::pair<int, GLuint>> loadedTiles;
		std::vector<std::pair<int, GLuint>> loadedBigTiles;
		std::vector<std::pair<int, GLuint>> loadedSprites;
		int selectedTileIdx = 0;
		int selectedBigTileIdx = 0;
		int selectedSpriteIdx = 0;
		uint voxelValue;

		State* stateHead;
		State* stateTail;
		State* stateCurrent;

		bool enabled = false;
		bool undoEnabled = true;
		uint allocatedUndo = 0;

		std::set<int3, std::less<>> editedBricks;
	};
}