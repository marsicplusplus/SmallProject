#pragma once
#include <future>

class FluidSimulator
{
public:
	FluidSimulator(World* world);
	~FluidSimulator();

	void Update(float deltaTime);
	void SetMaterialBlock(uint x, uint y, uint z, uint width, uint height, uint depth, uint amount, bool clear);
private:
	void InitCAPE(World* world, uint updateRate);
	void UpdateCAPE(float deltaTime);
	void CAPEThread(float deltaTime);

	bool capeRunning = false;
	std::future<void> capeThread;
	CAPE* cape;
};

