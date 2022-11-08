#pragma once

namespace Tmpl8
{

	class MyGame : public Game
	{

		void PrintDebug();
		void PrintStats();
	public:
		// game flow methods
		void Init();
		void Tick(float deltaTime);
		void Shutdown() { /* implement if you want to do something on exit */ }
		// input handling 
		void MouseUp(int button) { /* implement if you want to do something on mouse up */ }
		void MouseDown(int button) { /* implement if you want to do something on mouse down */ }
		void MouseMove(int x, int y) { /* implement if you want to do something on mouse move */  }
		void KeyUp(int key) { /* implement if you want to do something on key up */ }
		void KeyDown(int key) { /* implement if you want to do something on key down */  }
		// data members
		void DumpScreenBuffer();

		void HandleControls(float deltaTime);
		void PreRender();

		//wrapper for void(MyGame, int) fn
		static void IntArgFunction(std::function<void(MyGame&, int)> fn, MyGame& g, std::string s, int defaultarg);
	};

} // namespace Tmpl8