#pragma once

namespace Tmpl8
{
	class WorldEditor
	{
		typedef struct Selected
		{
			aabb box;
			float3 N;
		};

		// Taking inspiration from Goxel's gesture class
		enum GestureButton {
			GESTURE_LMB = 0,
			GESTURE_RMB
		};

		enum GestureState {
			GESTURE_POSSIBLE = 0,
			GESTURE_START,
			GESTURE_UPDATE,
			GESTURE_END
		};

		typedef struct Gesture
		{
			GestureButton button;
			GestureState state;
		};

	public:
		WorldEditor();
		void MouseMove(int x, int y);
		void MouseDown(int button);
		void MouseUp(int button);

	private:
		void UpdateSelectedBrick();
		void RemoveBrick();
		void AddBrick();
		bool HitWorldGrid(const float3 O, const float3 D);

		Selected selectedBricks;
		int2 mousePos;
		Gesture gesture;
		PAYLOAD* tempBricks = 0;
		uint* tempGrid;

	};
}

