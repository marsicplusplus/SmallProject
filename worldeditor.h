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

	public:
		void MouseMove(int x, int y);
		void MouseDown(int button);

	private:
		void UpdateSelectedVoxels();
		bool HitWorldGrid(const float3 O, const float3 D);

		Selected selectedVoxels;
		int2 mousePos;

	};
}

