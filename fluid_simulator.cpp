#include "precomp.h"
#include "CAPE.h"
#include "fluid_simulator.h"

FluidSimulator::FluidSimulator(World* world)
{
	InitCAPE(world, 100);
}

FluidSimulator::~FluidSimulator()
{
	delete cape;
}

void FluidSimulator::Update(float deltaTime)
{
	cape->Tick(deltaTime);
}

void FluidSimulator::SetMaterialBlock(uint x, uint y, uint z, uint width, uint height, uint depth, uint amount, bool clear)
{
	cape->SetMaterialBlock(x, y, z, width, height, depth, amount, clear);
}

void FluidSimulator::InitCAPE(World* world, uint updateRate)
{
	cape = new CAPE();
	cape->Initialise(world, updateRate);
	cout << "Initialised CAPE" << endl;
}

void FluidSimulator::UpdateCAPE(float deltaTime)
{
	//start new sim if done with previous update (no currently active)
	if(!capeThread.valid())
		capeThread = std::async(&FluidSimulator::CAPEThread, this, deltaTime);
	if (capeThread.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
		capeThread = std::async(&FluidSimulator::CAPEThread, this, deltaTime);
}

void FluidSimulator::CAPEThread(float deltaTime)
{
	cape->Tick(deltaTime);
}
