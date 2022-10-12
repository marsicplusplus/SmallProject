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
		enum class GestureButton {
			GESTURE_DEFAULT = -1,
			GESTURE_LMB = 0,
			GESTURE_RMB,
		};

		enum class GestureKey {
			GESTURE_DEFAULT = -1,
			GESTURE_CTRL = 0
		};

		enum class GestureState {
			GESTURE_POSSIBLE = 0,
			GESTURE_START,
			GESTURE_UPDATE,
			GESTURE_END
		};

		enum class GestureMode {
			GESTURE_DEFAULT = -1,
			GESTURE_ADD = 0,
			GESTURE_SUBTRACT
		};

		typedef struct Gesture
		{
			GestureButton button = GestureButton::GESTURE_DEFAULT;
			GestureKey key = GestureKey::GESTURE_DEFAULT;
			GestureState state = GestureState::GESTURE_POSSIBLE;
			GestureMode mode = GestureMode::GESTURE_ADD;
		};

	public:
		WorldEditor();
		void MouseMove(int x, int y);
		void MouseDown(int button);
		void MouseUp(int button);
		void KeyUp(int key);
		void KeyDown(int key);

	private:
		void UpdateSelectedBrick();
		void RemoveBrick();
		void AddBrick();
		bool HitWorldGrid(const float3 O, const float3 D);

		Selected selectedBricks;
		int2 mousePos;
		Gesture gesture;
		PAYLOAD* tempBricks = 0;
		uint* tempGrid = 0;
	};
}

