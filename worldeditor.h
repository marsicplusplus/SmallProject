#pragma once

namespace Tmpl8
{
	class WorldEditor
	{
		typedef struct Selected
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
			GESTURE_START,
			GESTURE_UPDATE
		};

		enum GestureMode {
			GESTURE_ADD = 0,
			GESTURE_REMOVE = 1 << 0,
			GESTURE_MULTI = 1 << 1
		};

		typedef struct Gesture
		{
			uint buttons = GestureButton::GESTURE_NO_BUTTONS;
			uint keys = GestureKey::GESTURE_NO_KEYS;
			uint state = GestureState::GESTURE_POSSIBLE;
			uint mode = GestureMode::GESTURE_ADD;
		};

	public:
		WorldEditor();
		void MouseMove(int x, int y);
		void MouseDown(int button);
		void MouseUp(int button);
		void KeyUp(int key);
		void KeyDown(int key);
		bool IsEnabled() { return enabled; }
		void Enable() { enabled = true; }
		void Disable() { enabled = false; }

	private:
		void UpdateSelectedBrick();
		void RemoveBrick();
		void AddBrick();
		void MultiAddRemove();
		void UpdateGestureMode();

		Selected selectedBricks;
		int2 mousePos;
		Gesture gesture;
		uint selectedButtons = GestureButton::GESTURE_NO_BUTTONS;
		uint selectedKeys = GestureKey::GESTURE_NO_KEYS;

		PAYLOAD* tempBricks = 0;
		uint* tempGrid = 0;
		std::vector<int> loadedTiles;
		int tileIdx;
		bool enabled = false;
	};
}

