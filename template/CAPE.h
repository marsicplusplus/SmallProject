#pragma once

namespace Tmpl8
{

	// CAPE - Cellular Automata Physics Engine
#define USECONCURRENCY 1
#define PRESSURE_ITERATIONS 8
#define DEBUG_MODE 1
#define CAPE_BRICKDIM 4
#define CAPE_GRIDWIDTH (MAPWIDTH / CAPE_BRICKDIM + 2)
#define CAPE_GRIDHEIGHT (MAPHEIGHT / CAPE_BRICKDIM + 2)
#define CAPE_GRIDDEPTH (MAPDEPTH / CAPE_BRICKDIM + 2)
#define CAPE_GRIDSIZE (CAPE_GRIDWIDTH * CAPE_GRIDHEIGHT * CAPE_GRIDDEPTH)
#define CAPE_BRICKSIZE (CAPE_BRICKDIM * CAPE_BRICKDIM * CAPE_BRICKDIM)
#define BIX(x, y, z) ((x)+(y) * (CAPE_GRIDWIDTH) + (z) * (CAPE_GRIDHEIGHT) * (CAPE_GRIDDEPTH))
#define MIN_BRICKMASS 0.001f
#define VELOCITY_DAMPENING 0.1f 
#define AL 0.5f //Advection limiter, prevents oscillations
#define INVAL (1.0f/AL) 


#define EVAPORATION 0.001f //amount of mass to remove per cell per second (helps clean up low mass cells that are not visible.)
#define MINRENDERMASS 0.001f //Dont render voxels below this much mass
#define GRAVITYENABLED 1
#define CELLSIZE 0.25f //In meters, essentially multiplier for gravity

//0.333 is maximum 100% save speed, use up to 1.0f for less clamping and faster flowing water (less viscous), but requires that an appropriate
//timestep is selected by the user, or simulation may blow up if local velocity becomes too high and negative mass is created
#define MAXV 1.0f //0.33f

	class CAPE
	{
	public:
		CAPE() {};
		~CAPE();

		World* world;

		//Run CA physics for this frame
		void Tick(float deltaTime);
		void Initialise(World* w, uint updateRate);
		void SetMaterialBlock(uint x, uint y, uint z, uint width, uint height, uint depth, float amount, bool clear = true);
		void ConvertToVoxels();
		void SetColorForCell(uint x, uint y, uint z, float timeStep);
		void AddMaterial(uint x, uint y, uint z, float amount);
		void ClearMaterial(uint x, uint y, uint z);

	private:
		//Configuration
		uint solveIterations = 1;
		float timeStep = 0.01;
		float brickTime = 0;
		float advectTime = 0;
		float velupdateTime = 0;
		float divergenceupdatetime = 0;
		float pressuresolvetime = 0;
		float pressuregradienttime = 0;
		float prevtime = 0;
		int updates = 0;

		Timer timer;
		float simulationTime = 0;
		int cellUpdates = 0;

		vector<float3> ga;

		uint GetBrickIDX(const uint x, const uint y, const uint z);
		float GetData(const uint x, const uint y, const uint z, vector<float*>& data);
		void SetData(const uint x, const uint y, const uint z, float v, vector<float*>& data);
		void AddData(const uint x, const uint y, const uint z, float v, vector<float*>& data);
		void EraseBrick(uint i);
		void FreeBrick(uint i);
		uint NewBrick(uint bx, uint by, uint bz);
		void CheckCompressMemory();
		void UpdateBricks();

		//Brick addresses, or uint max if brick empty
		uint* grid;
		uint bricks_alive = 0;
		uint bricks_allocated = 0;

		//total mass of active bricks
		vector<float> brick_m;
		vector<bool> brick_static;

		vector<uint> trash;
		vector<uint> brick_x;
		vector<uint> brick_y;
		vector<uint> brick_z;

		vector<float*> m_bricks; //Bricks containing all material data
		vector<float*> m0_bricks;
		vector<float*> p_bricks;
		vector<float*> p0_bricks;
		vector<float*> div_bricks;
		vector<float*> vx_bricks;
		vector<float*> vy_bricks;
		vector<float*> vz_bricks;
		vector<float*> vx0_bricks;
		vector<float*> vy0_bricks;
		vector<float*> vz0_bricks;

		float TotalDivergence();
		void PrintState();
		float MaterialChange(uint x, uint y, uint z, float vxl, float vxr, float vyl, float vyr, float vzl, float vzr);
		float IncomingMomentumX(uint x, uint y, uint z);
		float IncomingMomentumY(uint x, uint y, uint z);
		float IncomingMomentumZ(uint x, uint y, uint z);
		bool IsCellStatic(uint x, uint y, uint z);
		void CellVelocityUpdate(uint x, uint y, uint z, float timeStep);
		void SolvePressure(uint x, uint y, uint z, float timeStep);
		void CellDivergenceUpdate(uint x, uint y, uint z, float timeStep);
		float TotalMass();
		void MaterialAdvection(uint x, uint y, uint z, float timeStep);
		void PressureGradient(uint x, uint y, uint z, float timeStep);
		void RunOverAllBricks(CAPE* cape, void(CAPE::* func)(uint, uint, uint, float), float timeStep);
	};

}

