#pragma once
#include <CAPE.h>
#include <future>


#define THREADSAFEWORLD 1
#define SQR(x) ((x)*(x))
#define TILESIZE	8
#define TILESIZE2	(TILESIZE * TILESIZE)

#define OUTOFRANGE -99999

namespace Tmpl8
{
// Sprite system overview:
// The world contains a set of 0 or more sprites, typically loaded from .vox files.
// Sprites act like classic homecomputer sprites: they do not affect the world in
// any way, they just get displayed. Internally, this works as follows:
// 1. Before Game::Tick is executed:
//    - the world gets rendered by the GPU
//    - the sprites are then removed from the world
// 2. Game::Tick is now executed on a world that does not contain the sprites.
// 3. After Game::Tick completes:
//    - each sprite makes a backup of the voxels it overlaps
//    - the sprites are added back to the world
//    - the world is synchronized with the GPU for rendering in step 1.
// Consequence of this system is that even stationary sprites take time to process.

class SpriteFrame
{
	// fast sprite system: B&H'21
public:
	~SpriteFrame() { _aligned_free( buffer ); }
	PAYLOAD* buffer = 0;				// full frame buffer (width * height * depth)
	int3 size = make_int3( 0 );			// size of the sprite over x, y and z
	uint* drawPos = 0;					// compacted list of opaque sprite voxel positions
	PAYLOAD* drawVal = 0;				// compacted list of opaque sprite voxel values
	uint drawListSize = 0;				// number of voxels in buffer2
};

class Sprite
{
public:
	vector<SpriteFrame*> frame;			// sprite frames
	SpriteFrame* backup = 0;			// backup of pixels that the sprite overwrote
	int3 lastPos = make_int3( OUTOFRANGE );	// location where the backup will be restored to
	int3 currPos = make_int3( OUTOFRANGE );	// location where the sprite will be drawn
	mat4 transform = mat4::Identity();	// only 3x3 part is used; scale + rotation
	int currFrame = 0;					// frame to draw
	bool hasShadow = false;				// set to true to enable a drop shadow
	uint4* preShadow = 0;				// room for backup of voxels overwritten by shadow
	uint shadowVoxels = 0;				// size of backup voxel array
};

class SpriteManager
{
	// this class replaces the original World::sprite vector; it exists purely so that
	// we can create sprites at the global scope safely, should we so desire.
public:
	SpriteManager() = default;
	static SpriteManager* GetSpriteManager()
	{
		if (!spriteManager) spriteManager = new SpriteManager();
		return spriteManager;
	}
	uint LoadSprite( const char* voxFile, bool largeModel = false );
	void SaveSprite( const uint idx, const char* vxFile );
	uint CloneSprite( const uint idx );
	vector<Sprite*> sprite;				// list of loaded sprites
private:
	static inline SpriteManager* spriteManager = 0;
};

class Particles
{
public:
	Particles( int N )					// constructor
	{
		count = N, voxel = new uint4[N], backup = new uint4[N];
		memset( voxel, 0, N * sizeof( uint4 ) );
		memset( backup, 0, N * sizeof( uint4 ) );
		voxel[0].x = OUTOFRANGE;		// inactive by default
	}
	uint4* voxel = 0;					// particle positions & color
	uint4* backup = 0;					// backup of voxels overlapped by particles
	uint count = 0;						// particle count for the set
};

class ParticlesManager
{
	// this class replaces the original World::particles vector;
	// see SpriteManager for details.
public:
	ParticlesManager() = default;
	static ParticlesManager* GetParticlesManager()
	{
		if (!particlesManager) particlesManager = new ParticlesManager();
		return particlesManager;
	}
	uint CreateParticles( const uint count );
	// data members
	vector<Particles*> particles;		// list of particle sets
private:
	static inline ParticlesManager* particlesManager = 0;
};

// Tile system overview:
// The top-level grid / brick layout of the world (see below) fits well with the 
// classic concept of tiled graphics. A tile is simply an 8x8x8 or 16x16x16 chunk of 
// voxel data, which can be placed in the world at locations that are a multiple of 8 
// over x, y and z. Drawing a tile will thus simply overwrite the contents of a brick 
// (or 8 bricks, when using the larger 16x16x16 tiles).

class Tile
{
public:
	Tile() = default;
	Tile( const char* voxFile );
	PAYLOAD voxels[BRICKSIZE];			// tile voxel data
	uint zeroes;						// number of transparent voxels in the tile
};

class BigTile
{
public:
	BigTile() = default;
	BigTile( const char* voxFile );
	Tile tile[8];						// a big tile is just 2x2x2 tiles stored together
};

class TileManager
{
	// this class replaces the original World::tile and World::bigTile vectors;
	// see SpriteManager for details.
public:
	TileManager() = default;
	static TileManager* GetTileManager()
	{
		if (!tileManager) tileManager = new TileManager();
		return tileManager;
	}
	uint LoadTile( const char* voxFile );
	uint LoadBigTile( const char* voxFile );
	// data members
	vector<Tile*> tile;					// list of loaded tiles
	vector<BigTile*> bigTile;			// list of loaded big tiles
private:
	static inline TileManager* tileManager = 0;
};

// Voxel world data structure:
// The world consists of a 128x128x128 top-level grid. Each cell in this grid can
// either store a solid color, or the index of an 8x8x8 brick. Filling all cells with 
// brick indices yields the maximum world resolution of 1024x1024x01024.
// Voxels are 8-bit values. '0' is an empty voxel; all other colors are opaque. Voxel
// colors are 3-3-2 rgb values. Note that black voxels do not exist in this scheme.
// The data structure is mirrored to the GPU, with a delay of 1 frame (i.e., the GPU
// always renders the previous frame). 
// Furthermore, since the full dataset is a bit over 1GB, only changes are synced.
// the CPU to GPU communication consists of 8MB for the 128x128x128 top-level ints,
// plus up to 8192 changed bricks. If more changes are made per frame, these will
// be postponed to the next frame.

class World
{
public:
	//thread capeThread;
	// constructor / destructor
	CAPE* cape;
	World(const uint targetID);
	~World();
	// initialization
	void Clear();
	void Fill(const uint c);
	void DummyWorld();
	void LoadSky(const char* filename, const float scale = 1.0f);
	float3 SampleSky(const float3& D);
	void UpdateSkylights(); // updates the six skylight colors
	void ForceSyncAllBricks();
	void OptimizeBricks();
	// camera
	void SetCameraMatrix(const mat4& m) { camMat = m; }
	float3 GetCameraViewDir() { return make_float3(camMat[2], camMat[6], camMat[10]); }
	void SetCameraViewDir(const float3 D) { camMat[2] = D.x, camMat[6] = D.y, camMat[10] = D.z; }
	float3 GetCameraPos() { return make_float3(camMat[3], camMat[7], camMat[11]); }
	void SetCameraPos(const float3 P) { camMat[3] = P.x, camMat[7] = P.y, camMat[11] = P.z; }
	mat4& GetCameraMatrix() { return camMat; }
	float4* GetSkyLight() { return skyLight; }
	RenderParams& GetRenderParams() { return params; }
	DebugInfo& GetDebugInfo() { return debugInfo; }
	cl_event& GetRenderDoneEventHandle() { return renderDone; }
	cl_event& GetCommitDoneEventHandle() { return commitDone; }
	cl_event& GetCopyDoneEventHandle() { return copyDone; }
	Buffer* GetFrameBuffer() { return tmpFrame; }
	Buffer* GetAccumulatorBuffer() { return accumulator; }
	Buffer* GetLightsBuffer() { return lightsBuffer; }
	Buffer* GetDebugBuffer() { return debugBuffer; }
	void SetIsViewer(bool v) { viewer = v; }
	bool IsViewer() { return viewer; }
	Buffer** GetReservoirsBuffer() { return reservoirBuffers; }
	void SetViewerParamBuffer(Buffer* s) { viewerParamsBuffer = s; }
	void SetViewerPixelBuffer(Buffer* s) { viewerPixelBuffer = s; }
	Buffer* GetScreenBuffer() { return screen; }
	void SetScreenBuffer(Buffer* s) { screen = s; }
	uint GetTextureId() { return targetTextureID; }
	uint* GetGrid() { return grid; }
	PAYLOAD* GetBrick() { return brick; }
	uint* GetTrash() { return trash; }
	uint* GetZeroes() { return zeroes; }
	void SetLightsBuffer(Buffer* buffer) { lightsBuffer = buffer; };
	void SetReservoirBuffer(Buffer* buffer, int index) { reservoirBuffers[index] = buffer; }
	void Commit();
	void Render();
	float GetAverage(float* values, unsigned int numValues);
	float GetRenderTime() { return GetAverage(renderTimes, NumFrametimeSamples); }
	float GetAlbedoTime() { return GetAverage(albedoTimes, NumFrametimeSamples); }
	float GetCandidateTime() { return GetAverage(candidateTimes, NumFrametimeSamples); }
	float GetSpatialTime() { return GetAverage(spatialTimes, NumFrametimeSamples); }
	float GetShadingTime() { return GetAverage(shadingTimes, NumFrametimeSamples); }
	float GetFinalizingTime() { return GetAverage(finalizingTimes, NumFrametimeSamples);  }
	// high-level voxel access
	void Sphere( const float x, const float y, const float z, const float r, const uint c );
	void HDisc( const float x, const float y, const float z, const float r, const uint c );
	void Print( const char* text, const uint x, const uint y, const uint z, const uint c );
	void PrintZ( const char* text, const uint x, const uint y, const uint z, const uint c );
	uint CreateSprite( const int3 pos, const int3 size, const int frames );
	uint SpriteFrameCount( const uint idx );
	void MoveSpriteTo( const uint idx, const uint x, const uint y, const uint z );
	void TransformSprite( const uint idx, mat4 transform );
	void RemoveSprite( const uint idx );
	void StampSpriteTo( const uint idx, const uint x, const uint y, const uint z );
	void SetSpriteFrame( const uint idx, const uint frame );
	bool SpriteHit( const uint A, const uint B );
	void SetParticle( const uint set, const uint idx, const uint3 pos, const uint v );
	void DrawTile( const uint idx, const uint x, const uint y, const uint z );
	void DrawTiles( const char* tileString, const uint x, const uint y, const uint z );
	void DrawBigTile( const uint idx, const uint x, const uint y, const uint z );
	void DrawBigTiles( const char* tileString, const uint x, const uint y, const uint z );
	// inline ray tracing / cpu-only ray tracing / inline ray batch rendering
	uint TraceRay( float4 A, const float4 B, float& dist, float3& N, int steps );
	uint TraceBrick( float4 A, const float4 B, float& dist, float3& N, int steps );
	uint TraceTile(float4 A, const float4 B, const PAYLOAD* tile, float& dist, float3& N, int steps);
	uint TraceRay(float4 A, const float4 B, float& dist, float3& N, int steps, const PAYLOAD* oldBricks, const uint* oldGrid);
	void TraceRayToVoid( float4 A, const float4 B, float& dist, float3& N );
	Ray* GetBatchBuffer();
	Intersection* TraceBatch( const uint batchSize );
	Intersection* TraceBatchToVoid( const uint batchSize );
	// block scrolling
	void ScrollX( const int offset );
	void ScrollY( const int offset );
	void ScrollZ( const int offset );
	aabb& GetBounds() { return bounds; }

	// Lights
	void UpdateLights(float deltaTime);
	void FindLightsInWord(vector<Light> &ls);
	void SetupLightBuffer() { vector<Light> empty; SetupLightBuffer(empty); }
	void SetupLightBuffer(const vector<Light>& ls, int pos = 0);

	void AddLight(const int3 pos, const uint size, const uint c);

	void AddRandomLights(int numberOfLights);
	void RemoveRandomLights(int numberOfLights);
	void RemoveLight(const int3 pos);

	void SetUpMovingLights(int numberOfLights);
	void MoveLights();
	void PopLights(float deltaTime);
	uint MovingLightCount() { return movinglights.size(); }

	void InitReSTIR();
	float GetCurrentFrametime();

	vector<Tile*>& GetTileList() { return TileManager::GetTileManager()->tile; }
	vector<BigTile*>& GetBigTileList() { return TileManager::GetTileManager()->bigTile; }
	// convenient access to 'guaranteed to be instantiated' sprite, particle, tile lists
	vector<Sprite*>& GetSpriteList() { return SpriteManager::GetSpriteManager()->sprite; }
	vector<Particles*>& GetParticlesList() { return ParticlesManager::GetParticlesManager()->particles; }

private:
	// internal methods
	void EraseSprite( const uint idx );
	void DrawSprite( const uint idx );
	void DrawSpriteShadow( const uint idx );
	void RemoveSpriteShadow( const uint idx );
	void EraseParticles( const uint set );
	void DrawParticles( const uint set );
	void DrawTileVoxels( const uint cellIdx, const PAYLOAD* voxels, const uint zeroes );
	void SetupReservoirBuffers();

public:
	// low-level voxel access
	__forceinline PAYLOAD Get( const uint x, const uint y, const uint z)
	{
		// calculate brick location in top-level grid
		const uint bx = (x / BRICKDIM) & (GRIDWIDTH - 1);
		const uint by = (y / BRICKDIM) & (GRIDHEIGHT - 1);
		const uint bz = (z / BRICKDIM) & (GRIDDEPTH - 1);
		const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
		const uint cellValue = grid[cellIdx];
		if (IsSolidGridCell(cellValue) /* this is currently a 'solid' grid cell */) return cellValue >> 1;
		// calculate the position of the voxel inside the brick
		const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
		const uint brickBufferOffset = cellValue >> 1;
		return brick[brickBufferOffset * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM];
	}

	__forceinline void AddBrick(const uint bx, const uint by, const uint bz, const uint v /* actually an 8-bit value */)
	{
		if (bx >= GRIDWIDTH || by >= GRIDHEIGHT || bz >= GRIDDEPTH) return;
		const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
		uint cellValue = grid[cellIdx], brickBufferOffset = cellValue >> 1;

		// If we're adding a brick to a non-solid cell,
		// we need to free the brick and ensure that any lights in the cell
		// are removed from the light buffer
		if (!IsSolidGridCell(cellValue))
		{
			FreeBrick(brickBufferOffset);
			for (uint ly = 0; ly < BRICKDIM; ly++)
				for (uint lz = 0; lz < BRICKDIM; lz++)
					for (uint lx = 0; lx < BRICKDIM; lx++)
					{
						uint _v = Get(lx + bx * BRICKDIM, ly + by * BRICKDIM, lz + bz * BRICKDIM);

						// If we're replacing what was a light, remove it from the light buffer
						if (IsEmitter(_v))
						{
							RemoveLight(int3(lx + bx * BRICKDIM, ly + by * BRICKDIM, lz + bz * BRICKDIM));
						}
					}
		}
		else
		{
			// If we're replacing what was a brick light, remove it from the light buffer
			uint oldV = cellValue >> 1;
			if (IsEmitter(oldV))
			{
				RemoveLight(int3(bx * BRICKDIM, by * BRICKDIM, bz * BRICKDIM));
			}

		}

		// If we're adding a brick light, add it to the light buffer
		if (IsEmitter(v))
		{
			AddLight(int3(bx * BRICKDIM, by * BRICKDIM, bz * BRICKDIM), BRICKDIM, v);
		}

		Mark(0); // Mark to ensure new grid gets sent to GPU
		grid[cellIdx] = v << 1;
	}

	__forceinline void RemoveBrick(const uint bx, const uint by, const uint bz)
	{
		if (bx >= GRIDWIDTH || by >= GRIDHEIGHT || bz >= GRIDDEPTH) return;
		const uint cellIndex = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
		// obtain current brick identifier from top-level grid
		uint cellValue = grid[cellIndex], brickBufferOffset = cellValue >> 1;
		grid[cellIndex] = 0;	// brick just became completely zeroed; recycle

		if (!IsSolidGridCell(cellValue)) // If not solid/empty, free brick
		{
			zeroes[brickBufferOffset] = BRICKSIZE;
			FreeBrick(brickBufferOffset);
			for (uint ly = 0; ly < BRICKDIM; ly++)
				for (uint lz = 0; lz < BRICKDIM; lz++)
					for (uint lx = 0; lx < BRICKDIM; lx++)
					{
						uint v = Get(lx + bx * BRICKDIM, ly + by * BRICKDIM, lz + bz * BRICKDIM);

						// If we're removing individual voxels that were lights, remove them from the light buffer
						if (IsEmitter(v))
						{
							RemoveLight(int3(lx + bx * BRICKDIM, ly + by * BRICKDIM, lz + bz * BRICKDIM));
						}
					}
		}
		else
		{
			// If we're replacing what was a brick light, remove it from the light buffer
			uint oldV = cellValue >> 1;
			if (IsEmitter(oldV))
			{
				RemoveLight(int3(bx * BRICKDIM, by * BRICKDIM, bz * BRICKDIM));
			}
		}

		Mark(0); // Mark to ensure new grid gets sent to GPU
	}

	__forceinline int SplitSolidBrick(uint brickColor, uint cellIndex)
	{
		const uint newBrickBufferOffset = NewBrick();

		#if BRICKDIM == 8 && PAYLOADSIZE == 1
			// fully unrolled loop for writing the 512 bytes needed for a single brick, faster than memset
			const __m256i brickColor8 = _mm256_set1_epi8( static_cast<char>(brickColor) );
			__m256i* d8 = (__m256i*)(brick + newBrickBufferOffset * BRICKSIZE);
			d8[0] = brickColor8, d8[1] = brickColor8, d8[2] = brickColor8, d8[3] = brickColor8;
			d8[4] = brickColor8, d8[5] = brickColor8, d8[6] = brickColor8, d8[7] = brickColor8;
			d8[8] = brickColor8, d8[9] = brickColor8, d8[10] = brickColor8, d8[11] = brickColor8;
			d8[12] = brickColor8, d8[13] = brickColor8, d8[14] = brickColor8, d8[15] = brickColor8;
		#elif BRICKDIM == 8 && PAYLOADSIZE == 2
			// fully unrolled loop for writing 1KB needed for a single brick, faster than memset
			const __m256i brickColor16 = _mm256_set1_epi16( static_cast<short>(brickColor) );
			__m256i* d = (__m256i*)(brick + newBrickBufferOffset * BRICKSIZE);
			d[0] = brickColor16, d[1] = brickColor16, d[2] = brickColor16, d[3] = brickColor16;
			d[4] = brickColor16, d[5] = brickColor16, d[6] = brickColor16, d[7] = brickColor16;
			d[8] = brickColor16, d[9] = brickColor16, d[10] = brickColor16, d[11] = brickColor16;
			d[12] = brickColor16, d[13] = brickColor16, d[14] = brickColor16, d[15] = brickColor16;
			d[16] = brickColor16, d[17] = brickColor16, d[18] = brickColor16, d[19] = brickColor16;
			d[20] = brickColor16, d[21] = brickColor16, d[22] = brickColor16, d[23] = brickColor16;
			d[24] = brickColor16, d[25] = brickColor16, d[26] = brickColor16, d[27] = brickColor16;
			d[28] = brickColor16, d[29] = brickColor16, d[30] = brickColor16, d[31] = brickColor16;
		#elif BRICKDIM == 8 && PAYLOADSIZE == 4
			// fully unrolled loop for writing 2KB needed for a single brick, faster than memset
			const __m256i brickColor32 = _mm256_set1_epi32( static_cast<uint>(brickColor) );
			__m256i* d = (__m256i*)(brick + newBrickBufferOffset * BRICKSIZE);
			d[0] = brickColor32, d[1] = brickColor32, d[2] = brickColor32, d[3] = brickColor32;
			d[4] = brickColor32, d[5] = brickColor32, d[6] = brickColor32, d[7] = brickColor32;
			d[8] = brickColor32, d[9] = brickColor32, d[10] = brickColor32, d[11] = brickColor32;
			d[12] = brickColor32, d[13] = brickColor32, d[14] = brickColor32, d[15] = brickColor32;
			d[16] = brickColor32, d[17] = brickColor32, d[18] = brickColor32, d[19] = brickColor32;
			d[20] = brickColor32, d[21] = brickColor32, d[22] = brickColor32, d[23] = brickColor32;
			d[24] = brickColor32, d[25] = brickColor32, d[26] = brickColor32, d[27] = brickColor32;
			d[28] = brickColor32, d[29] = brickColor32, d[30] = brickColor32, d[31] = brickColor32;
			d[32] = brickColor32, d[33] = brickColor32, d[34] = brickColor32, d[35] = brickColor32;
			d[36] = brickColor32, d[37] = brickColor32, d[38] = brickColor32, d[39] = brickColor32;
			d[40] = brickColor32, d[41] = brickColor32, d[42] = brickColor32, d[43] = brickColor32;
			d[44] = brickColor32, d[45] = brickColor32, d[46] = brickColor32, d[47] = brickColor32;
			d[48] = brickColor32, d[49] = brickColor32, d[50] = brickColor32, d[51] = brickColor32;
			d[52] = brickColor32, d[53] = brickColor32, d[54] = brickColor32, d[55] = brickColor32;
			d[56] = brickColor32, d[57] = brickColor32, d[58] = brickColor32, d[59] = brickColor32;
			d[60] = brickColor32, d[61] = brickColor32, d[62] = brickColor32, d[63] = brickColor32;
		#else
			// TODO: Generic case
		#endif

		zeroes[newBrickBufferOffset] = (brickColor == 0) * BRICKSIZE;
		// Update the grid, as the brick inside the grid is now an offset into the brick buffer
		// rather than a solid color
		grid[cellIndex] = (newBrickBufferOffset << 1) | 1;
		return newBrickBufferOffset;
	}

	__forceinline void RemoveBigTile(const uint bx, const uint by, const uint bz)
	{
		for (int x = 0; x < 2; x++)
			for (int y = 0; y < 2; y++)
				for (int z = 0; z < 2; z++)
					RemoveBrick(bx + x, by + y, bz + z);
	}


	__forceinline void Set( const uint x, const uint y, const uint z, const uint v /* actually an 8-bit value */)
	{
		// calculate brick location in top-level grid
		uint bx = x / BRICKDIM;
		uint by = y / BRICKDIM;
		uint bz = z / BRICKDIM;

		if (bx >= GRIDWIDTH || by >= GRIDHEIGHT || bz >= GRIDDEPTH)
			return;		//Way to prevent this branching?

		const uint cellIndex = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
		uint cellValue = grid[cellIndex];
		uint brickColor = cellValue >> 1;
		uint brickBufferOffset = cellValue >> 1;

		if (IsSolidGridCell(cellValue))
		{
			// No change to the actual brick color, return
			if (v == brickColor)
			{
				return;
			}
			brickBufferOffset = SplitSolidBrick(brickColor, cellIndex);
			cellValue = grid[cellIndex];
		}


		// calculate the position of the voxel inside the brick
		const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
		const uint voxelIdx = brickBufferOffset * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM; //Precalculate this?
		const uint originalVoxel = brick[voxelIdx];

		// If we're replacing what was a light, remove it from the light buffer
		if (IsEmitter(originalVoxel))
		{
			RemoveLight(int3(x, y, z));
		}

		// If we're adding a light, add it to the light buffer
		if (IsEmitter(v))
		{
			AddLight(int3(x, y, z), 1, v);
		}


		int zeroChange = (originalVoxel != 0 && v == 0) - (originalVoxel == 0 && v != 0);
		zeroes[brickBufferOffset] += zeroChange;
		if (zeroes[brickBufferOffset] < BRICKSIZE)
		{
			brick[voxelIdx] = v; 
			Mark(brickBufferOffset);
		}
		else
		{
			grid[cellIndex] = 0;
			FreeBrick(brickBufferOffset);
		}
	}

	__forceinline void SetFromSprite(const uint x, const uint y, const uint z, const uint v /* actually an 8-bit value */)
	{
		// calculate brick location in top-level grid
		uint bx = x / BRICKDIM;
		uint by = y / BRICKDIM;
		uint bz = z / BRICKDIM;

		if (bx >= GRIDWIDTH || by >= GRIDHEIGHT || bz >= GRIDDEPTH)
			return;		//Way to prevent this branching?

		const uint cellIndex = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
		uint cellValue = grid[cellIndex];
		uint brickColor = cellValue >> 1;
		uint brickBufferOffset = cellValue >> 1;

		if (IsSolidGridCell(cellValue))
		{
			// No change to the actual brick color, return
			if (v == brickColor)
			{
				return;
			}
			brickBufferOffset = SplitSolidBrick(brickColor, cellIndex);
			cellValue = grid[cellIndex];
		}

		// calculate the position of the voxel inside the brick
		const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
		const uint voxelIdx = brickBufferOffset * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM; //Precalculate this?
		const uint originalVoxel = brick[voxelIdx];

		int zeroChange = (originalVoxel != 0 && v == 0) - (originalVoxel == 0 && v != 0);
		zeroes[brickBufferOffset] += zeroChange;
		if (zeroes[brickBufferOffset] < BRICKSIZE)
		{
			brick[voxelIdx] = v;
			Mark(brickBufferOffset);
		}
		else
		{
			grid[cellIndex] = 0;
			FreeBrick(brickBufferOffset);
		}
	}

	//Temp Getter to allow easy access to world data from CAPE
	Buffer* GetBrickBuffer() { return brickBuffer; }
	Buffer* GetZeroesBuffer() { return zeroesBuffer; }
	cl_mem GetGridMap() { return gridMap; }

	void Mark(const uint idx)
	{
	#if THREADSAFEWORLD
		// be careful, setting a bit in an array is not thread-safe without _interlockedbittestandset
		_interlockedbittestandset((LONG*)modified + (idx >> 5), idx & 31);
	#else
		modified[idx >> 5] |= 1 << (idx & 31);
	#endif
	}

	void UnMark(const uint idx)
	{
	#if THREADSAFEWORLD
		// be careful, resetting a bit in an array is not thread-safe without _interlockedbittestandreset
		_interlockedbittestandreset((LONG*)modified + (idx >> 5), idx & 31);
	#else
		modified[idx >> 5] &= 0xffffffffu - (1 << (idx & 31));
	#endif
	}
public:
	uint NewBrick()
	{
	#if THREADSAFEWORLD
		// get a fresh brick from the circular list in a thread-safe manner and without false sharing
		const uint trashItem = InterlockedAdd( &trashTail, 31 ) - 31;
		return trash[trashItem & (BRICKCOUNT - 1)];
	#else
		// slightly faster to not prevent false sharing if we're doing single core updates only
		return trash[trashTail++ & (BRICKCOUNT - 1)];
	#endif
	}
	void FreeBrick( const uint idx )
	{
	#if THREADSAFEWORLD
		// thread-safe access of the circular list
		const uint trashItem = InterlockedAdd( &trashHead, 31 ) - 31;
		trash[trashItem & (BRICKCOUNT - 1)] = idx;
	#else
		// for single-threaded code, a stepsize of 1 maximizes cache coherence.
		trash[trashHead++ & (BRICKCOUNT - 1)] = idx;
	#endif
	}

private:


	bool IsDirty( const uint idx ) { return (modified[idx >> 5] & (1 << (idx & 31))) > 0; }
	bool IsDirty32( const uint idx ) { return modified[idx] != 0; }
	void ClearMarks32( const uint idx ) { modified[idx] = 0; }
	void ClearMarks() { memset( modified, 0, (BRICKCOUNT / 32) * 4 ); }
	// helpers
	__forceinline static void StreamCopy( __m256i* dst, const __m256i* src, const uint bytes )
	{
		// https://stackoverflow.com/questions/2963898/faster-alternative-to-memcpy
		assert( (bytes & 31) == 0 );
		if (!CPUCaps::HW_AVX2)
		{
			// fallback: no AVX2, use SSE 4.2
			uint N = bytes / 16;
			const __m128i* src4 = (__m128i*)src;
			__m128i* dst4 = (__m128i*)dst;
			for (; N > 0; N--, src4++, dst4++)
			{
				const __m128i d = _mm_stream_load_si128( src4 );
				_mm_stream_si128( dst4, d );
			}
		}
		else
		{
			// AVX2 path - this version: CO'21
			uint N = bytes / 32;
			constexpr uint registers = 8;
			uint unalignedStep = N % registers;
			for (; N > 0 && unalignedStep > 0; N--, unalignedStep--, src++, dst++)
			{
				const __m256i d = _mm256_stream_load_si256( src );
				_mm256_stream_si256( dst, d );
			}
			static_assert(registers == 8);
			for (; N > 0; N -= registers)
			{
				// Based on https://stackoverflow.com/questions/62419256/how-can-i-determine-how-many-avx-registers-my-processor-has
				const __m256i d0 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d0 );
				const __m256i d1 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d1 );
				const __m256i d2 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d2 );
				const __m256i d3 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d3 );
				const __m256i d4 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d4 );
				const __m256i d5 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d5 );
				const __m256i d6 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d6 );
				const __m256i d7 = _mm256_stream_load_si256( src++ );
				_mm256_stream_si256( dst++, d7 );
			}
		}
	}
	void StreamCopyMT( __m256i* dst, __m256i* src, const uint bytes );
	// helper class for multithreaded memcpy
	class CopyJob : public Job
	{
	public:
		void Main() { World::StreamCopy( dst, src, N * 32 ); }
		__m256i* dst, * src;
		uint N;
	};
	// data members
	aabb bounds;
	mat4 camMat;						// camera matrix to be used for rendering
	uint* grid = 0, * gridOrig = 0;		// pointer to host-side copy of the top-level grid
#if ONEBRICKBUFFER == 1
	Buffer* brickBuffer;				// OpenCL buffer for the bricks
#else
	Buffer* brickBuffer[4];				// OpenCL buffers for the bricks
#endif
	PAYLOAD* brick = 0;					// pointer to host-side copy of the bricks
	uint* modified = 0;					// bitfield to mark bricks for synchronization
	volatile inline static LONG trashHead = BRICKCOUNT;	// thrash circular buffer tail
	volatile inline static LONG trashTail = 0;	// thrash circular buffer tail
	bool viewer = false;
	uint* trash = 0;					// indices of recycled bricks
	uint* zeroes = 0;
	Buffer* zeroesBuffer;
	Buffer* screen = 0;					// OpenCL buffer that encapsulates the target OpenGL texture
	uint targetTextureID = 0;			// OpenGL render target
	int prevFrameIdx = 0;				// index of the previous frame buffer that will be used for TAA
	Buffer* debugBuffer = 0;
	Buffer* paramBuffer = 0;			// OpenCL buffer that stores renderer parameters
	Buffer* historyTAA[2] = { 0 };			// OpenCL buffers for history data (previous frame)
	Buffer* tmpFrame = 0;				// OpenCL buffer to store rendered frame in linear color space
	Buffer* accumulator = 0;
	Buffer* sky = 0;					// OpenCL buffer for a HDR skydome
	Buffer* blueNoise = 0;				// blue noise data
	Buffer* lightsBuffer = 0;
	Buffer* reservoirBuffers[2] = { 0 };
	Buffer* primaryHitBuffer[2] = { 0 };
	Buffer* viewerParamsBuffer;
	Buffer* viewerPixelBuffer;
	int2 skySize;						// size of the skydome bitmap
	RenderParams params;				// CPU-side copy of the renderer parameters
	DebugInfo debugInfo;
	Kernel* albedoRender;
	Kernel* perPixelLightSampling;
	Kernel* spatialResampling;
	Kernel* finalizeSimple;
	Kernel* currentRenderer;
	Kernel* committer;					// render kernel and commit kernel
	Kernel* finalizer, * unsharpen;		// TAA finalization kernels
	Kernel* accumulatorFinalizer;		// Path tracer finalization kernels
	Kernel* rendererTAA, * rendererNoTAA, * rendererRIS, * rendererPathTracer;
	Kernel* batchTracer;				// ray batch tracing kernel for inline tracing
	Kernel* batchToVoidTracer;			// ray batch tracing kernel for inline tracing from solid to void
#if MORTONBRICKS == 1
	Kernel* encodeBricks;				// reorganizes brick data on the GPU
#endif
	cl_mem uberGrid = 0;				// 32x32x32 top-level grid (device-side only)
	Kernel* uberGridUpdater;			// build a 32x32x32 top-level grid over the brickmap
	cl_event ubergridDone;				// for profiling
	cl_event copyDone, commitDone;		// events for queue synchronization
	public:
	cl_event renderDone;				// event used for profiling
	cl_event albedoRenderDone;
	cl_event candidateAndTemporalResamplingDone;
	cl_event spatialResamplingDone;
	cl_event shadingDone;
	cl_event finalizingDone;
	private:

	constexpr static uint NumFrametimeSamples = 100;
	float renderTimes[NumFrametimeSamples];					// render time for the previous frame (in seconds)
	float albedoTimes[NumFrametimeSamples];
	float candidateTimes[NumFrametimeSamples];
	float spatialTimes[NumFrametimeSamples];
	float shadingTimes[NumFrametimeSamples];
	float finalizingTimes[NumFrametimeSamples];
	uint64_t currentFrame = 0;
	uint tasks = 0;						// number of changed bricks, to be passed to commit kernel
	bool copyInFlight = false;			// flag for skipping async copy on first iteration
	bool commitInFlight = false;		// flag to make next commit wait for previous to complete
	cl_mem devmem = 0;					// device-side commit buffer
	cl_mem gridMap;						// device-side 3D image or buffer for top-level
	Surface* font;						// bitmap font for print command
	bool firstFrame = true;				// for doing things in the first frame
	float4 skyLight[6];					// integrated light for the 6 possible normals
	unordered_map<uint, uint> defaultVoxel;
	unordered_map<uint, uint4> movinglights;
	public:
	bool lightsAreMoving = false;
	bool poppingLights = false;
};

} // namespace Tmpl8
