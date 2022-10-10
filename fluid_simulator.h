#pragma once
#include <future>

class FluidSimulator
{
public:
	FluidSimulator(World* world);
	~FluidSimulator();

	void Update(float deltaTime);
	void SetMaterialBlock(const float x, const float y, const float z, const float width, const float height, const float depth, const float amount, const bool clear);
private:
	void InitCAPE(World* world, uint updateRate);
	void UpdateCAPE(float deltaTime);
	void CAPEThread(float deltaTime);

	bool capeRunning = false;
	std::future<void> capeThread;
	CAPE* cape;
};

