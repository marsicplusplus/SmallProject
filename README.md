# Voxel ReSTIR
Thesis INFOMGMT'21 Xander Hermans: The effectiveness of the ReSTIR technique when ray tracing a voxel world. [paper](https://github.com/xanderhermans/WrldTmpl8/blob/github/INFOMGMT_21_Thesis_Xander_Hermans.pdf)

![img](https://raw.githubusercontent.com/xanderhermans/WrldTmpl8/github/screenbufferdata/headerpic.png)

## Abstract
Recent work by Bitterli et al. combines previous techniques such as
weighted reservoir sampling and resampled importance sampling
into an algorithm called ReSTIR capable of rendering of scenes
illuminated by many dynamic lights in real-time. We introduce
Voxel ReSTIR, an algorithm derived from ReSTIR tailored to voxel
worlds, and compare it to ReSTIR. We explore the performance of
the algorithms in different settings. We show that that both our
point light and area light algorithms perform better in terms of
speed and quality than ReSTIR in the same setting. Our algorithm
is suitable for GPU implementation and we are able to achieve near
real-time performance on consumer level hardware.

## ReSTIR
The ReSTIR algorithm by Bitterli et al. uses Resampled Importance Sampling (RIS) to sample from a pool of emitters.
The algorithm incorporates Weighted Reservoir Sampling (WRS) to resample a number of candidates per pixel with only constant memory per pixel required.
After resampling the candidates the chosen candidate is tested for visibility given the shading point. 
If the candidate is not visible from the shading point the reservoir is discarded in which case the weight and length of the reservoir are set to 0.
The algorithm then applies spatial and temporal resampling.
The temporal resampling combines the reservoir for the current pixel with the reservoir from the pixel at the same location in last frame. 
When combining a reservoir with that of a previous frame the weight can grow unbounded.
To combat this the weight of a reservoir from a previous frame is clamped to 20 times the weight of the reservoir after the initial candidate sampling.
The spatial resampling picks 5 neighbouring pixels from a 30 pixel radius disc around the current pixel and combines the reservoirs from those neighbours with the reservoir of the current pixel.
By employing spatialtemporal resampling we effectively re-use the information gained from candidate sampling for many more candidates than just the candidates sampled for the current pixel.
This approach suffers from bias as the spatialtemporal resampling can lead to invalid candidates for a given pixel due to difference in visibility or geometry.

There are two ways in which ReSTIR deals with this bias.
1. Accept and manage the bias by using a heuristic to ignore reservoirs during spatialtemporal resampling.
Re-use visibility and assume each previous or neighbour reservoir candidate is visible from the shading point.
The result has a visible darkening bias but does not require extra rays to be traced.
2. Combine reservoirs using Multiple Importance Sampling (MIS) weights and include visibility in this equation.
The result does not suffer from the darkening bias but requires an additional ray per neighbour reservoir used in spatial resampling to be traced.

We used the ReSTIR implementation [GfxExp](https://github.com/shocker-0x15/GfxExp) by shocker-0x15 to compare our implementation to.

![restir](https://github.com/xanderhermans/WrldTmpl8/raw/github/restir_letters1.png)

## Voxel ReSTIR
Voxel ReSTIR is very similar to ReSTIR as described above. 
We decided to create Voxel ReSTIR based on the biased variant of the ReSTIR algorithm.
We chose the biased variant since our target is games and therefore real-time framerates.

The Voxel ReSTIR algorithm shares the same stages (candidate sampling, temporal resampling, spatial resampling, shading).
The difference between Voxel ReSTIR and ReSTIR is the way emitting primitives are sampled. 
Where normal ReSTIR uses triangle emitters Voxel ReSTIR uses voxel emitters.
This means that every emitter always has a positive Solid Angle (SA) greater than 0 to a shading point in the scene.
This property is exploited when picking a point on the emitting voxel candidate to sample from, as we employ importance sampling to choose a side of the voxel to sample from based on the solid angle.
It also increases the quality of our candidate sampling because every candidate has a positive non-zero potential contribution.
Our target Probability Density Function (PDF) used when candidate sampling is the potential contribution of the emitter much like the target PDF used by ReSTIR.
The difference is that since we importance sample the sides of the voxel, this probability is included in the potential contribution.

The temporal and spatial resampling stages in Voxel ReSTIR are similar to the stages in ReSTIR. 
ReSTIR uses a depth and geometry heuristic to determine if the reservoir for a given former or neighbouring pixel should be ignored.
The heuristic is as follows: For a shading point p and a former or neighbouring pixel q, accept q iff the depth difference with relation to the camera of p and q is no greater than 10% and the angle between p and q is no greater than 25 degrees. (distance(depth(p), depth(q)) <= 0.1) && (dot(normal(p), normal(q)) >= 0.5).
For Voxel ReSTIR we decided to keep the heuristic similar. While this meant no change for the depth heuristic the geometry heuristic can be simplified to normal(p) == normal(q) or side_of_voxel(p) == side_of_voxel(q) since we only have 6 distinct normal vectors per voxel and each forms an angle of at least 90 degrees with another. 

The Voxel ReSTIR shading stage has no ReSTIR specific differences.

The voxel emitters can be considered as area lights as described above but they can also be considered as point lights in our implementation.
Doing so achieves a considerable speedup of the candidate sampling stage which is beneficial for reaching real-time framerates.
Sadly, using point lights comes at the cost of hard shades and simpler highlights.

![pointlights](https://github.com/xanderhermans/WrldTmpl8/raw/github/point_lights_mountain2.png)

Our Voxel ReSTIR algoritm can deal with voxel emitters of arbitrary sizes.
This lead to experimenting with grouping of homogeneous spaces of emitting voxels into a single large emitting cube.
In our implementation we only group such spaces if they happen to be aligned as brick of 8x8x8 voxels in our grid since the logic for these bricks was already used in WrldTmpl8.
By grouping emitters into a single emitter we reduce the number of candidates which in turn increases the quality of the candidate sampling as there are fewer candidates to choose from, which in terms means that the candidate sampling more likely to closely represent the pool of candidates.

Due to the way the colour information per voxel is stored in WrldTmpl8 we only have 8 bits for emission strength of an emitter.
Currently this means an emitter can have strength [0..15] but this could easily be changed to an arbiratry range with 8 bit precision or 8 bit index into a collection of materials.

The skydome sampling in this implementation is separate from ReSTIR and is not considered during candidate sampling but rather always sampled during the shading stage.

![bricks_voxels](https://github.com/xanderhermans/WrldTmpl8/raw/github/screenbufferdata/results/bricks_voxels.png)
A comparison between voxel emitters grouped as bricks of 8x8x8 emitting voxels and individual voxel emitters.

![restir_voxel_restir_reference](https://github.com/xanderhermans/WrldTmpl8/raw/github/screenbufferdata/results/gfxexp_sa_gt/edited/gfxexp_sa_gt/flyingapartments1.png)
A comparison between ReSTIR, Voxel ReSTIR and the ground truth using the same settings and scene.

![point_area_reference](https://github.com/xanderhermans/WrldTmpl8/raw/github/screenbufferdata/results/point_sa_gt/edited/point_sa_gt/flyingapartments2.png)
A comparison between considering voxel emitters as point lights or considering voxel emitters as area lights.

![sa32_sa16_sa8_sa1](https://github.com/xanderhermans/WrldTmpl8/raw/github/screenbufferdata/results/sa_32_16_8_1/edited/sa_32_16_8_1/mountain1.png)
Varying number of candidates using Voxel ReSTIR with voxel emitters as area lights.

For more details and images refer to the [paper](https://github.com/xanderhermans/WrldTmpl8/blob/github/INFOMGMT_21_Thesis_Xander_Hermans.pdf).

## Implementation details
The Voxel ReSTIR project is implemented in the MyGame subproject. Make sure to set this project as startup project when running from Visual Studio.

At the top of [common.h](https://github.com/xanderhermans/WrldTmpl8/blob/github/template/common.h) are the definitions which control the behaviour of the implementation.
This is where the skydome sampling can be turned On or Off as well as various ReSTIR specific parameters.

This repository comes with 4 scenes in [scene_export](https://github.com/xanderhermans/WrldTmpl8/tree/github/scene_export).
The logic for generating, loading and saving these scenes is in [mygame.cpp:MyGame::Init()](https://github.com/xanderhermans/WrldTmpl8/blob/github/mygame.cpp).
This is also where the voxel world is scanned for emitting voxels which are added to the pool of lights.

The light manager in [light_manager.h](https://github.com/xanderhermans/WrldTmpl8/blob/github/light_manager.h) can be used to control the lights (emitting voxels) in the scene.
This class also has logic to move lights around without readding them to the pool of lights such that ReSTIR can take advantage of the temporal information. 
Note that this leaves visible artifacts unless both spatial and temporal resampling are both enabled.

The controls are defined in [mygame.cpp:MyGame::HandleControls(float deltaTime)](https://github.com/xanderhermans/WrldTmpl8/blob/github/mygame.cpp). \
W A S D for panning the camera around. \
Arrow keys for rotating the camera. \
R and F for panning up and down. \
Q for dumping the current screenbuffer as raster of float4 values to "screenbufferdata/tmpl8/". \
SPACE to toggle accumulating pixel information over time to converge when stationary. \
C and V to toggle spatial and temporal resampling respectively.

The OpenCL implementation of Voxel ReSTIR can be found in [ris.cl](https://github.com/xanderhermans/WrldTmpl8/blob/github/cl/ris.cl) which depends on [phat.cl](https://github.com/xanderhermans/WrldTmpl8/blob/github/cl/phat.cl), [reservoir.cl](https://github.com/xanderhermans/WrldTmpl8/blob/github/cl/reservoir.cl), [spatialSampling.cl](https://github.com/xanderhermans/WrldTmpl8/blob/github/cl/spatialSampling.cl) and [tools.cl](https://github.com/xanderhermans/WrldTmpl8/blob/github/cl/tools.cl).

The project has been maintained using Visual Studio 17.1.1 on Windows 10 (21H1) using an AMD r5 5600x cpu and Nvidia GTX 1080 gpu. 

## References
[RIS] (https://scholarsarchive.byu.edu/etd/663) \
[WRS] (https://doi.org/10.1016/j.ipl.2005.11.003) \
[ReSTIR] (https://doi.org/10.1145/3386569.3392481) \
[WrldTmpl8] (https://github.com/jbikker/WrldTmpl8) \
[GfxExp] (https://github.com/shocker-0x15/GfxExp)

Xander Hermans 2022