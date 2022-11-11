#pragma once
#include <fluid_simulator.h>

namespace Tmpl8
{

	struct PathPoint { float3 O, D; };

	class Lighthouse : public Game
	{
	public:
		Lighthouse();
		void InitialiseLighthouseScenario();
		void SetStaticBlock(uint x0, uint y0, uint z0, uint w, uint h, uint d, uint v);

		// game flow methods
		void Init();
		void HandleInput(float deltaTime);
		void Tick(float deltaTime);
		void Shutdown()
		{
			FILE* f = fopen("camera.dat", "wb"); // save camera
			fwrite(&D, 1, sizeof(D), f);
			fwrite(&O, 1, sizeof(O), f);
			fclose(f);
		}
		// input handling
		void MouseUp(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseDown(int button) { /* implement if you want to detect mouse button presses */ }
		void MouseMove(int x, int y) { mousePos.x = x, mousePos.y = y; }
		void KeyUp(int key) { /* implement if you want to handle keys */ }
		void KeyDown(int key) { /* implement if you want to handle keys */ }
		float3 CatmullRom(const float3& p0, const float3& p1, const float3& p2, const float3& p3)
		{
			const float3 c = 2 * p0 - 5 * p1 + 4 * p2 - p3, d = 3 * (p1 - p2) + p3 - p0;
			return 0.5f * (2 * p1 + ((p2 - p0) * t) + (c * t * t) + (d * t * t * t));
		}
		// data members
		int2 mousePos;
		float3 ballPos = make_float3(300, 100, 300);
		float3 ballVel = make_float3(0.3f, 0, 0.5f);
		// spline path data
		vector<PathPoint> splinePath;
		int pathPt = 1; // spline path vertex
		float t = 0; // spline path t (0..1)
		// camera
		float3 D = make_float3(0.931, -0.1081, -0.348);
		float3 O = make_float3(298, 530, 596);

		int lighthouseSprite, frame = 0;
		bool keyPressed[0xFF];
	};

} // namespace Tmpl8