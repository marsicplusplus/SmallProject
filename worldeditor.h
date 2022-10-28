#pragma once

namespace Tmpl8
{
	class WorldEditor
	{

		typedef struct State;
		struct State {
			PAYLOAD* newBricks;
			PAYLOAD* oldBricks;

			uint* newGridVals;
			uint* oldGridVals;

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
			GESTURE_TILE = 1,
			GESTURE_BIG_TILE = 2
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
		void Enable() { UpdateSelectedBrick(); enabled = true; }
		void Disable() { ResetEditor(); enabled = false; }
		void ResetEditor();
		void RenderGUI();

	private:
		void UpdateSelectedBrick();
		void AddBrick(int bx, int by, int bz);
		void RemoveBrick(int bx, int by, int bz);
		void MultiAddRemove();
		void UpdateGestureMode();
		void CheckMemoryAllowance();
		void Redo();
		void Undo();
		void SaveState();
		void DeleteState(State* state);
		bool CreateNewState();
		void LoadTiles();
		void UpdateEditedBricks(int bx, int by, int bz);

		// Input and Gesture 
		int2 mousePos;
		Gesture gesture;
		uint selectedButtons = GestureButton::GESTURE_NO_BUTTONS;
		uint selectedKeys = GestureKey::GESTURE_NO_KEYS;
		Selected selectedBricks;

		// Temporary buffers to hold previous state
		PAYLOAD* tempBricks = 0;
		uint* tempGrid = 0;
		BrickInfo* tempBrickInfo = 0;

		std::vector<std::pair<int, GLuint>> loadedTiles;
		std::vector<std::pair<int, GLuint>> loadedBigTiles;
		int selectedTileIdx = 0;
		int selectedBigTileIdx = 0;

		State* stateHead;
		State* stateTail;
		State* stateCurrent;
		bool enabled = false;
		bool undoEnabled = true;
		uint allocatedUndo = 0;

		std::vector<int3> editedBricks;
	};
}

