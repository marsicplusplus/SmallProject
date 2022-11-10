#include "precomp.h"
#include "lighthouse.h"

Game* CreateGame() { return new Lighthouse(); }

Tmpl8::Lighthouse::Lighthouse() {}

void Lighthouse::SetStaticBlock(uint x0, uint y0, uint z0, uint w, uint h, uint d, uint v)
{
	if (x0 + w > MAPWIDTH || y0 + h > MAPWIDTH || z0 + d > MAPWIDTH)
	{
		cout << "Block outside of grid" << endl;
		return;
	}
	for (uint x = x0; x < x0 + w; x++)
		for (uint y = y0; y < y0 + h; y++)
			for (uint z = z0; z < z0 + d; z++)
				Plot(x, y, z, v);
}

void Lighthouse::InitialiseLighthouseScenario()
{
	lighthouseSprite = LoadSprite("assets/lighthouse.vox");
	MoveSpriteTo(lighthouseSprite, make_int3(450, 500, 500));

	World& world = *GetWorld();
	for (int lx = 18; lx < 25; ++lx) {
		for (int ly = 62; ly < 71; ++ly) {
			world.Set(450 + lx, 500 + ly, 500 + 60, WHITE | (255 << 16));
		}
	}
	for (int lx = 18; lx < 25; ++lx) {
		for (int ly = 62; ly < 71; ++ly) {
			world.Set(450 + lx, 500 + ly, 500 + 68, WHITE | (255 << 16));
		}
	}
	for (int ly = 62; ly < 70; ++ly) {
		for (int lz = 61; lz < 68; ++lz) {
			world.Set(450 + 25, 500 + ly, 500 + lz, WHITE | (255 << 16));
		}
	}
	for (int ly = 62; ly < 70; ++ly) {
		for (int lz = 61; lz < 68; ++lz) {
			world.Set(450 + 17, 500 + ly, 500 + lz, WHITE | (255 << 16));
		}
	}
	SetStaticBlock(450 + 114, 500 + 5, 500 + 8, 2, 3, 2, WHITE | 100 << 16);
	SetStaticBlock(450 + 88, 500 + 5, 500 + 26, 2, 3, 2, WHITE | 100 << 16);
	SetStaticBlock(450 + 112, 500 + 5, 500 + 67, 2, 3, 2, WHITE | 100 << 16);
	SetStaticBlock(450 + 87, 500 + 5, 500 + 105, 2, 3, 2, WHITE | 100 << 16);
	SetStaticBlock(400, 400, 400, 300, 103, 300, 0X00f | 13 << 12);
}

static bool shouldDumpBuffer = false;
static bool takeScreenshot = false;
static bool useSpatialResampling = USESPATIAL;
static bool useTemporalResampling = USETEMPORAL;
static bool skyDomeSampling = true;


// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Lighthouse::Init()
{
	World& world = *GetWorld();
	ShowCursor(false);

	InitialiseLighthouseScenario();

	/* Overwrite defaults for ReSTIR */
	RenderParams& params = world.GetRenderParams();
	params.numberOfLights = 0;
	params.accumulate = false;
	params.editorEnabled = false;
	params.spatial = useSpatialResampling;
	params.temporal = useTemporalResampling;
	params.spatialTaps = SPATIALTAPS;
	params.spatialRadius = SPATIALRADIUS;
	params.numberOfCandidates = NUMBEROFCANDIDATES;
	params.numberOfMaxTemporalImportance = TEMPORALMAXIMPORTANCE;
	params.skyDomeSampling = skyDomeSampling;
	world.GetDebugInfo().counter = 0;

	world.OptimizeBricks(); //important to recognize bricks
	vector<Light> vls;
	world.FindLightsInWord(vls);
	world.SetupLightBuffer(vls);
	skyDomeLightScale = 0.0f;
	skyDomeImage = "assets/sky_21.hdr";
}

// -----------------------------------------------------------
// HandleInput: reusable input handling / free camera
// -----------------------------------------------------------
void Lighthouse::HandleInput(float deltaTime)
{
	// free cam controls
	float3 tmp(0, 1, 0), right = normalize(cross(tmp, D)), up = cross(D, right);
	float speed = deltaTime * 0.1f;
	if (GetAsyncKeyState('W')) O += speed * D; else if (GetAsyncKeyState('S')) O -= speed * D;
	if (GetAsyncKeyState('A')) O -= speed * right; else if (GetAsyncKeyState('D')) O += speed * right;
	if (GetAsyncKeyState('R')) O += speed * up; else if (GetAsyncKeyState('F')) O -= speed * up;
	if (GetAsyncKeyState(VK_LEFT)) D = normalize(D - right * 0.025f * speed);
	if (GetAsyncKeyState(VK_RIGHT)) D = normalize(D + right * 0.025f * speed);
	if (GetAsyncKeyState(VK_UP)) D = normalize(D - up * 0.025f * speed);
	if (GetAsyncKeyState(VK_DOWN)) D = normalize(D + up * 0.025f * speed);
	if (GetAsyncKeyState(VK_SPACE))
	{
		SetSpriteFrame(lighthouseSprite, (frame = frame + 1) % 5);
	}
	if (GetAsyncKeyState(VK_SHIFT) && !keyPressed[VK_SHIFT])
	{
		runCAPESimulation = !runCAPESimulation;
		keyPressed[VK_SHIFT] = true;
	}
	else if (!GetAsyncKeyState(VK_SHIFT))
	{
		keyPressed[VK_SHIFT] = false;
	}

	if (GetAsyncKeyState('G') && !keyPressed['G'])
	{
		World& world = *GetWorld();
		WorldEditor& worldEditor = *world.getWorldEditor();
		if (worldEditor.IsEnabled())
		{
			ShowCursor(false);
			worldEditor.Disable();
		}
		else
		{
			ShowCursor(true);
			worldEditor.Enable();
		}
		keyPressed['G'] = true;
	}
	else if (!GetAsyncKeyState('G'))
	{
		keyPressed['G'] = false;
	}


#if 1
	// enable to set spline path points using P key
	static bool pdown = false;
	static FILE* pf = 0;
	if (!GetAsyncKeyState('P')) pdown = false; else
	{
		if (!pdown) // save a point for the spline
		{
			if (!pf) pf = fopen("spline.bin", "wb");
			fwrite(&O, 1, sizeof(O), pf);
			float3 t = O + D;
			fwrite(&t, 1, sizeof(t), pf);
		}
		pdown = true;
	}
	LookAt(O, O + D);
#else
	// playback of recorded spline path
	const size_t N = splinePath.size();
	PathPoint p0 = splinePath[(pathPt + (N - 1)) % N], p1 = splinePath[pathPt];
	PathPoint p2 = splinePath[(pathPt + 1) % N], p3 = splinePath[(pathPt + 2) % N];
	LookAt(CatmullRom(p0.O, p1.O, p2.O, p3.O), CatmullRom(p0.D, p1.D, p2.D, p3.D));
	if ((t += deltaTime * 0.0005f) > 1) t -= 1, pathPt = (pathPt + 1) % N;
#endif
}

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Lighthouse::Tick(float deltaTime)
{
	// update camera
	HandleInput(deltaTime);
}
