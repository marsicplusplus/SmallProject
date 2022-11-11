#define DEBUGCOLOR (float3)(1.0, 0.0, 1.0)
#define DEBUGCOLOR2 (float3)(0.0, 1.0, 1.0)
#define DBX 500
#define DBY 500

float SphericalTheta(const float3 v)
{
	return acos(clamp(v.z, -1.f, 1.f));
}

float SphericalPhi(const float3 v)
{
	const float p = atan2(v.y, v.x);
	return (p < 0) ? (p + 2 * PI) : p;
}

uint WangHash(uint s) { s = (s ^ 61) ^ (s >> 16), s *= 9, s = s ^ (s >> 4), s *= 0x27d4eb2d, s = s ^ (s >> 15); return s; }
uint RandomInt(uint* s) { *s ^= *s << 13, * s ^= *s >> 17, * s ^= *s << 5; return *s; }
float RandomFloat(uint* s) { return RandomInt(s) * 2.3283064365387e-10f; }

float3 DiffuseReflectionCosWeighted(const float r0, const float r1, const float3 N)
{
	const float3 T = normalize(cross(N, fabs(N.y) > 0.99f ? (float3)(1, 0, 0) : (float3)(0, 1, 0)));
	const float3 B = cross(T, N);
	const float term1 = TWOPI * r0, term2 = sqrt(1 - r1);
	float c, s = sincos(term1, &c);
	return (c * term2 * T) + (s * term2) * B + sqrt(r1) * N;
}

// https://www.scratchapixel.com/code.php?id=34&origin=/lessons/3d-basic-rendering/global-illumination-path-tracing
float3 UniformSampleHemisphere(const float r1, const float r2, const float3 N)
{
	float3 Nt, Nb;
	if (fabs(N.x) > fabs(N.y))
		Nt = (float3)(N.z, 0, -N.x) / sqrt(N.x * N.x + N.z * N.z);
	else
		Nt = (float3)(0, -N.z, N.y) / sqrt(N.y * N.y + N.z * N.z);
	Nb = cross(N, Nt);


	// cos(theta) = r1 = y
	// cos^2(theta) + sin^2(theta) = 1 -> sin(theta) = srtf(1 - cos^2(theta))
	float sinTheta = sqrt(1 - r1 * r1);
	float phi = TWOPI * r2;
	float x = sinTheta * cos(phi);
	float z = sinTheta * sin(phi);
	return (float3)(
		x * Nb.x + r1 * N.x + z * Nt.x,
		x * Nb.y + r1 * N.y + z * Nt.y,
		x * Nb.z + r1 * N.z + z * Nt.z
		);
}

float blueNoiseSampler(const __global uint* blueNoise, int x, int y, int sampleIndex, int sampleDimension)
{
	// Adapated from E. Heitz. Arguments:
	// sampleIndex: 0..255
	// sampleDimension: 0..255
	x &= 127, y &= 127, sampleIndex &= 255, sampleDimension &= 255;
	// xor index based on optimized ranking
	int rankedSampleIndex = (sampleIndex ^ blueNoise[sampleDimension + (x + y * 128) * 8 + 65536 * 3]) & 255;
	// fetch value in sequence
	int value = blueNoise[sampleDimension + rankedSampleIndex * 256];
	// if the dimension is optimized, xor sequence value based on optimized scrambling
	value ^= blueNoise[(sampleDimension & 7) + (x + y * 128) * 8 + 65536];
	// convert to float and return
	float retVal = (0.5f + value) * (1.0f / 256.0f) /* + noiseShift (see LH2) */;
	if (retVal >= 1) retVal -= 1;
	return retVal;
}

// tc ∈ [-1,1]² | fov ∈ [0, π) | d ∈ [0,1] -  via https://www.shadertoy.com/view/tt3BRS
float3 PaniniProjection(float2 tc, const float fov, const float d)
{
	const float d2 = d * d;
	{
		const float fo = PI * 0.5f - fov * 0.5f;
		const float f = cos(fo) / sin(fo) * 2.0f;
		const float f2 = f * f;
		const float b = (native_sqrt(max(0.f, (d + d2) * (d + d2) * (f2 + f2 * f2))) - (d * f + f)) / (d2 + d2 * f2 - 1);
		tc *= b;
	}
	const float h = tc.x, v = tc.y, h2 = h * h;
	const float k = h2 / ((d + 1) * (d + 1)), k2 = k * k;
	const float discr = max(0.f, k2 * d2 - (k + 1) * (k * d2 - 1));
	const float cosPhi = (-k * d + native_sqrt(discr)) / (k + 1.f);
	const float S = (d + 1) / (d + cosPhi), tanTheta = v / S;
	float sinPhi = native_sqrt(max(0.f, 1 - cosPhi * cosPhi));
	if (tc.x < 0.0) sinPhi *= -1;
	const float s = native_rsqrt(1 + tanTheta * tanTheta);
	return (float3)(sinPhi, tanTheta, cosPhi) * s;
}

__constant float halton[32] = {
	0, 0, 0.5f, 0.333333f, 0, 0.6666666f, 0.75f, 0.111111111f, 0, 0.44444444f,
	0.5f, 0.7777777f, 0.25f, 0.222222222f, 0.75f, 0.55555555f, 0, 0.88888888f,
	0.5f, 0.03703703f, 0.25f, 0.37037037f, 0.75f, 0.70370370f, 0.125f, 0.148148148f,
	0.625f, 0.481481481f, 0.375f, 0.814814814f, 0.875f, 0.259259259f
};

// produce a camera ray direction for a position in screen space
float3 GenerateCameraRay(const float2 pixelPos, __constant struct RenderParams* params)
{
#if TAA
	const uint px = (uint)pixelPos.x;
	const uint py = (uint)pixelPos.y;
	const uint h = (params->frame + (px & 3) + 4 * (py & 3)) & 15;
	const float2 uv = (float2)(
		(pixelPos.x + halton[h * 2 + 0]) * params->oneOverRes.x,
		(pixelPos.y + halton[h * 2 + 1]) * params->oneOverRes.y
		);
#else
	const float2 uv = (float2)(
		pixelPos.x * params->oneOverRes.x,
		pixelPos.y * params->oneOverRes.y
		);
#endif
#if PANINI
	const float3 V = PaniniProjection((float2)(uv.x * 2 - 1, (uv.y * 2 - 1) * ((float)SCRHEIGHT / SCRWIDTH)), PI / 5, 0.15f);
	// multiply by improvised camera matrix
	return V.z * normalize((params->p1 + params->p2) * 0.5f - params->E) +
		V.x * normalize(params->p1 - params->p0) + V.y * normalize(params->p2 - params->p0);
#else
	const float3 P = params->p0 + (params->p1 - params->p0) * uv.x + (params->p2 - params->p0) * uv.y;
	return normalize(P - params->E);
#endif
}

// sample the HDR sky dome texture (bilinear)
float3 SampleSky(const float3 T, __global float4* sky, uint w, uint h)
{
	const float u = w * SphericalPhi(T) * INV2PI - 0.5f;
	const float v = h * SphericalTheta(T) * INVPI - 0.5f;
	const float fu = u - floor(u), fv = v - floor(v);
	const int iu = (int)u, iv = (int)v;
	const uint idx1 = (iu + iv * w) % (w * h);
	const uint idx2 = (iu + 1 + iv * w) % (w * h);
	const uint idx3 = (iu + (iv + 1) * w) % (w * h);
	const uint idx4 = (iu + 1 + (iv + 1) * w) % (w * h);
	const float4 s =
		sky[idx1] * (1 - fu) * (1 - fv) + sky[idx2] * fu * (1 - fv) +
		sky[idx3] * (1 - fu) * fv + sky[idx4] * fu * fv;
	return s.xyz;
}

bool HitAABB(const float3 O, const float3 D, const float3 boxMin, const float3 boxMax, bool useWinding, float* distance)
{
	float3 center = (boxMax + boxMin) / 2;
	float3 radius = (boxMax - boxMin) / 2;

	// Move to the box's reference frame. This is unavoidable and un-optimizable.
	float3 rayOrigin = O - center;


	// Winding direction: -1 if the ray starts inside of the box (i.e., and is leaving), +1 if it is starting outside of the box
	float3 res = fabs(rayOrigin) * (1.0 / radius);
	float max = res.x;
	if (res.y > max) max = res.y;
	if (res.z > max) max = res.z;
	float winding = (max < 1.0) ? -1.0 : 1.0;


	// We'll use the negated sign of the ray direction in several places, so precompute it.
	// The sign() instruction is fast...but surprisingly not so fast that storing the result
	// temporarily isn't an advantage.
	float3 sgn = -sign(D);

	// Ray-plane intersection. For each pair of planes, choose the one that is front-facing
	// to the ray and compute the distance to it.
	float3 distanceToPlane = radius * (useWinding ? winding : -1.0) * sgn - rayOrigin;
	distanceToPlane *= (1.0 / D);

	// Perform all three ray-box tests and cast to 0 or 1 on each axis. 
	// Use a macro to eliminate the redundant code (no efficiency boost from doing so, of course!)
	// Could be written with 
#   define TEST(U, VW)\
         /* Is there a hit on this axis in front of the origin? Use multiplication instead of && for a small speedup */\
         (distanceToPlane.U >= 0.0) && \
         /* Is that hit within the face of the box? */\
         all(isless(fabs(rayOrigin.VW + D.VW * distanceToPlane.U), radius.VW))

	int3 test = (int3)(TEST(x, yz), TEST(y, zx), TEST(z, xy));

	// CMOV chain that guarantees exactly one element of sgn is preserved and that the value has the right sign
	sgn = test.x ? (float3)(sgn.x, 0.0, 0.0) : (test.y ? (float3)(0.0, sgn.y, 0.0) : (float3)(0.0, 0.0, test.z ? sgn.z : 0.0));
#   undef TEST

	// At most one element of sgn is non-zero now. That element carries the negative sign of the 
	// ray direction as well. Notice that we were able to drop storage of the test vector from registers,
	// because it will never be used again.


	// Normal must face back along the ray. If you need
	// to know whether we're entering or leaving the box, 
	// then just look at the value of winding. If you need
	// texture coordinates, then use box.invDirection * hitPoint.
	float3 normal = sgn;

	// Mask the distance by the non-zero axis
	// Dot product is faster than this CMOV chain, but doesn't work when distanceToPlane contains nans or infs. 
	//
	float dist = (sgn.x != 0.0) ? distanceToPlane.x : ((sgn.y != 0.0) ? distanceToPlane.y : distanceToPlane.z);
	if (dist > *distance && *distance != 0)
		return false;
	//printf(" Distance: %f \n", *distance);

	*distance = dist;
	bool hit = (sgn.x != 0) || (sgn.y != 0) || (sgn.z != 0);
	return hit;
}

bool HitSelectedBox(const float3 O, const float3 D, const float3 bboxMin, const float3 bboxMax, const float wireBoxWidth, float* distance)
{
	float closeDistance = *distance;
	float3 wireMin = bboxMin;
	float3 wireMax = bboxMax + (float3)(1.0, 1.0, 1.0);

	if (HitAABB(O, D, wireMin, wireMax, true, &closeDistance))
	{
		const float3 shadingPoint = D * closeDistance + O;
		float delta = wireBoxWidth;
		bool hitx = fabs(shadingPoint.x - wireMin.x) < delta || fabs(wireMax.x - shadingPoint.x) < delta;
		bool hity = fabs(shadingPoint.y - wireMin.y) < delta || fabs(wireMax.y - shadingPoint.y) < delta;
		bool hitz = fabs(shadingPoint.z - wireMin.z) < delta || fabs(wireMax.z - shadingPoint.z) < delta;


		if (((hitx || hity) && hitz) ||
			((hity || hitz) && hitx) ||
			((hitx || hitz) && hity))
		{
			*distance = closeDistance;
			return true;
		}
	}

	float farDistance = *distance;
	if (HitAABB(O, D, wireMin, wireMax, false, &farDistance))
	{
		const float3 shadingPoint = D * farDistance + O;
		float delta = wireBoxWidth;
		bool hitx = fabs(shadingPoint.x - wireMin.x) < delta || fabs(wireMax.x - shadingPoint.x) < delta;
		bool hity = fabs(shadingPoint.y - wireMin.y) < delta || fabs(wireMax.y - shadingPoint.y) < delta;
		bool hitz = fabs(shadingPoint.z - wireMin.z) < delta || fabs(wireMax.z - shadingPoint.z) < delta;


		if (((hitx || hity) && hitz) ||
			((hity || hitz) && hitx) ||
			((hitx || hitz) && hity))
		{
			*distance = farDistance;
			return true;
		}
	}

	return false;
}

// Check if we've hit the world grid using the Ray-Box Intersection Algorithm from https://jcgt.org/published/0007/03/04/paper-lowres.pdf
bool HitWorldGrid(const float3 O, const float3 D)
{

	float3 worldMin = (0, 0, 0);
	float3 worldMax = (1024, 1024, 1024);
	float distance = 0;
	// This allows us to ensure we're always looking inside the world grid
	bool useWinding = false;
	if (!HitAABB(O, D, worldMin, worldMax, useWinding, &distance)) return false;

	//printf(/*"Shading Point X: %f, Shading Point Y: %f,  Shading Point Z: %f, */ "Distance: % f", /*shadingPoint.x, shadingPoint.y, shadingPoint.z,*/ distance);

	float3 hitpoint = O + D * distance;
	float delta = 0.1;

	// Check if we've hit a section of the grid
	//bool hitx = hitpoint.x - floor(hitpoint.x) < delta || ceil(hitpoint.x) - hitpoint.x < delta;
	//bool hity = hitpoint.y - floor(hitpoint.y) < delta || ceil(hitpoint.y) - hitpoint.y < delta;
	//bool hitz = hitpoint.z - floor(hitpoint.z) < delta || ceil(hitpoint.z) - hitpoint.z < delta;

	bool hitx = ((int)floor(hitpoint.x) % 8 == 0 || (int)ceil(hitpoint.x) % 8 == 0) /*&& (hitpoint.x - floor(hitpoint.x) < delta || ceil(hitpoint.x) - hitpoint.x < delta)*/;
	bool hity = ((int)floor(hitpoint.y) % 8 == 0 || (int)ceil(hitpoint.y) % 8 == 0) /*&& (hitpoint.y - floor(hitpoint.y) < delta || ceil(hitpoint.y) - hitpoint.y < delta)*/;
	bool hitz = ((int)floor(hitpoint.z) % 8 == 0 || (int)ceil(hitpoint.z) % 8 == 0) /*&& (hitpoint.z - floor(hitpoint.z) < delta || ceil(hitpoint.z) - hitpoint.z < delta)*/;

	if (((hitx || hity) && hitz) ||
		((hity || hitz) && hitx) ||
		((hitx || hitz) && hity))
	{
		return true;
	}

	return false;

}

uint GetAlpha(uint v)
{
	return (v & 0xf000) >> 12;
}

float GetAlphaf(uint v)
{
	uint alphaChannel = (v & 0xf000) >> 12;
	return alphaChannel / (float)0xF;
}

// convert a voxel color to floating point rgb
float3 ToFloatRGB(const uint v)
{
#if PAYLOADSIZE == 1
	return (float3)((v >> 5) * (1.0f / 7.0f), ((v >> 2) & 7) * (1.0f / 7.0f), (v & 3) * (1.0f / 3.0f));
#else
	return (float3)(((v >> 8) & 15) * (1.0f / 15.0f), ((v >> 4) & 15) * (1.0f / 15.0f), (v & 15) * (1.0f / 15.0f));
#endif
}

uint FromFloatRGBA(const float4 rgba, uint emitStrength)
{
	return (emitStrength << 16) +
		(convert_uint(rgba.a * 15) << 12) + 
		(convert_uint(rgba.r * 15) << 8) + 
		(convert_uint(rgba.g * 15) << 4) + 
		(convert_uint(rgba.b * 15));
}

// ACES filmic tonemapping, via https://www.shadertoy.com/view/3sfBWs
float3 ACESFilm(const float3 x)
{
	float a = 2.51f, b = 0.03f, c = 2.43f, d = 0.59f, e = 0.14f;
	return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0f, 1.0f);
}

// Reinhard2 tonemapping, via https://www.shadertoy.com/view/WdjSW3
float3 Reinhard2(const float3 x)
{
	const float3 L_white = (float3)4.0f;
	return (x * (1.0f + x / (L_white * L_white))) / (1.0f + x);
}

// https://twitter.com/jimhejl/status/633777619998130176
float3 ToneMapFilmic_Hejl2015(const float3 hdr, float whitePt)
{
	float4 vh = (float4)(hdr, whitePt);
	float4 va = 1.425f * vh + 0.05f;
	float4 vf = (vh * va + 0.004f) / (vh * (va + 0.55f) + 0.0491f) - 0.0821f;
	return vf.xyz / vf.www;
}

// Linear to SRGB, also via https://www.shadertoy.com/view/3sfBWs
float3 LessThan(const float3 f, float value)
{
	return (float3)(
		(f.x < value) ? 1.0f : 0.0f,
		(f.y < value) ? 1.0f : 0.0f,
		(f.z < value) ? 1.0f : 0.0f);
}
float3 LinearToSRGB(const float3 rgb)
{
#if 0
	return sqrt(rgb);
#else
	const float3 _rgb = clamp(rgb, 0.0f, 1.0f);
	return mix(pow(_rgb * 1.055f, (float3)(1.f / 2.4f)) - 0.055f,
		_rgb * 12.92f, LessThan(_rgb, 0.0031308f)
	);
#endif
}

// conversions between RGB and YCoCG for TAA
float3 RGBToYCoCg(const float3 RGB)
{
	const float3 rgb = min((float3)(4), RGB); // clamp helps AA for strong HDR
	const float Y = dot(rgb, (float3)(1, 2, 1)) * 0.25f;
	const float Co = dot(rgb, (float3)(2, 0, -2)) * 0.25f + (0.5f * 256.0f / 255.0f);
	const float Cg = dot(rgb, (float3)(-1, 2, -1)) * 0.25f + (0.5f * 256.0f / 255.0f);
	return (float3)(Y, Co, Cg);
}

float3 YCoCgToRGB(const float3 YCoCg)
{
	const float Y = YCoCg.x;
	const float Co = YCoCg.y - (0.5f * 256.0f / 255.0f);
	const float Cg = YCoCg.z - (0.5f * 256.0f / 255.0f);
	return (float3)(Y + Co - Cg, Y + Cg, Y - Co - Cg);
}

// sample a color buffer with bilinear interpolation
float3 bilerpSample(const __global float4* buffer, const float u, const float v)
{
	// do a bilerp fetch at u_prev, v_prev
	float fu = u - floor(u), fv = v - floor(v);
	int iu = (int)u, iv = (int)v;
	if (iu == 0 || iv == 0 || iu >= SCRWIDTH - 1 || iv >= SCRHEIGHT - 1) return (float3)0;
	int iu0 = clamp(iu, 0, SCRWIDTH - 1), iu1 = clamp(iu + 1, 0, SCRWIDTH - 1);
	int iv0 = clamp(iv, 0, SCRHEIGHT - 1), iv1 = clamp(iv + 1, 0, SCRHEIGHT - 1);
	float4 p0 = (1 - fu) * (1 - fv) * buffer[iu0 + iv0 * SCRWIDTH];
	float4 p1 = fu * (1 - fv) * buffer[iu1 + iv0 * SCRWIDTH];
	float4 p2 = (1 - fu) * fv * buffer[iu0 + iv1 * SCRWIDTH];
	float4 p3 = fu * fv * buffer[iu1 + iv1 * SCRWIDTH];
	return (p0 + p1 + p2 + p3).xyz;
}

// inefficient morton code for bricks (3 bit per axis, x / y / z)
int Morton3Bit(const int x, const int y, const int z)
{
	return (x & 1) + 2 * (y & 1) + 4 * (z & 1) +
		4 * (x & 2) + 8 * (y & 2) + 16 * (z & 2) +
		16 * (x & 4) + 32 * (y & 4) + 64 * (z & 4);
}
int Morton3(const int xyz)
{
	return (xyz & 1) + ((xyz & 8) >> 2) + ((xyz & 64) >> 4) +
		((xyz & 2) << 2) + (xyz & 16) + ((xyz & 128) >> 2) +
		((xyz & 4) << 4) + ((xyz & 32) << 2) + (xyz & 256);
}

// building a normal from an axis and a ray direction
float3 VoxelNormal(const uint side, const float3 D)
{
	if (side == 0) return (float3)(D.x > 0 ? -1 : 1, 0, 0);
	if (side == 1) return (float3)(0, D.y > 0 ? -1 : 1, 0);
	if (side == 2) return (float3)(0, 0, D.z > 0 ? -1 : 1);
}

int getRandomIndex(const int length, uint* seed)
{
	return clamp((int)(RandomFloat(seed) * length), 0, length - 1);
}

uint3 indexToCoordinates(const uint index)
{
	const uint y = index / (MAPWIDTH * MAPDEPTH);
	const uint z = (index / MAPWIDTH) % MAPDEPTH;
	const uint x = index % MAPDEPTH;
	return (uint3)(x, y, z);
}

uint Get(const uint x, const uint y, const uint z,
	__read_only image3d_t grid,
	__global const PAYLOAD* brick)
{
	// calculate brick location in top-level grid
	const uint bx = (x / BRICKDIM) & (GRIDWIDTH - 1);
	const uint by = (y / BRICKDIM) & (GRIDHEIGHT - 1);
	const uint bz = (z / BRICKDIM) & (GRIDDEPTH - 1);
	//const uint cellIdx = bx + bz * GRIDWIDTH + by * GRIDWIDTH * GRIDDEPTH;
	const uint g = read_imageui(grid, (int4)(bx, bz, by, 0)).x;
	if ((g & 1) == 0 /* this is currently a 'solid' grid cell */) return g >> 1;
	// calculate the position of the voxel inside the brick
	const uint lx = x & (BRICKDIM - 1), ly = y & (BRICKDIM - 1), lz = z & (BRICKDIM - 1);
	return brick[(g >> 1) * BRICKSIZE + lx + ly * BRICKDIM + lz * BRICKDIM * BRICKDIM];
}

uint3 GridCoordinatesFromHit(const int side, const float3 D, const float dist, const float3 O)
{
	const float3 N = VoxelNormal(side, D);
	const float3 hitPoint = D * dist + O;
	const float3 coordInsideVoxel = hitPoint - N * 0.1;
	return convert_uint3(coordInsideVoxel);
}

bool Uint3equals(const uint3 v1, const uint3 v2)
{
	return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

// v2.x is in [v1.x .. v1.x + size -1] and v2.y is in [v1.y .. v1.y + size -1] and v2.z is in [v1.z .. v1.z + size -1]
// simply put is v2 in the cube v1 to v1 + size.
bool Uint3CubeEquals(const uint3 v1, const uint3 v2, const uint size)
{
	return v1.x <= v2.x && v1.x + size > v2.x
		&& v1.y <= v2.y && v1.y + size > v2.y
		&& v1.z <= v2.z && v1.z + size > v2.z;
}

float3 ShadingPoint(const struct CLRay* ray, const float3 O)
{
	return ray->rayDirection * ray->distance + O;
}

int2 ReprojectWorldPoint(__constant struct RenderParams* params, const float3 pixelPos)
{
	const float dl = dot(pixelPos, params->Nleft.xyz) - params->Nleft.w;
	const float dr = dot(pixelPos, params->Nright.xyz) - params->Nright.w;
	const float dt = dot(pixelPos, params->Ntop.xyz) - params->Ntop.w;
	const float db = dot(pixelPos, params->Nbottom.xyz) - params->Nbottom.w;
	const float u_prev = SCRWIDTH * (dl / (dl + dr));
	const float v_prev = SCRHEIGHT * (dt / (dt + db));
	return (int2)(round(u_prev), round(v_prev));
}

// min inclusive, max exclusive
bool IsInBounds(const int2 coords, const int2 _min, const int2 _max)
{
	return coords.x >= _min.x && coords.y >= _min.y && coords.x < _max.x&& coords.y < _max.y;
}

//#define HALFNUMBEROFVOXELSIDES 3
#define NUMBEROFVOXELSIDES 6
// center of voxel(s), shading point, size of voxel(s) is the number of voxels in 1 dimension of the voxel cluster(Brick), ptr to seed, inverse probability for the point
float3 RandomPointOnVoxel(const float3 center, const float3 shadingPoint, int emitterDimensionSize, uint* seedptr, float* invProbability, float3* Nlight)
{
	const float size = convert_float(emitterDimensionSize);
	const float3 normals[NUMBEROFVOXELSIDES] =
	{
		(float3)(1.0, 0.0, 0.0),
		(float3)(0.0, 1.0, 0.0),
		(float3)(0.0, 0.0, 1.0),
		(float3)(-1.0, 0.0, 0.0),
		(float3)(0.0, -1.0, 0.0),
		(float3)(0.0, 0.0, -1.0),
	};
	const int t_indices[NUMBEROFVOXELSIDES] =
	{
		1, 2, 0, 4, 5, 3
	};
	const int b_indices[NUMBEROFVOXELSIDES] =
	{
		2, 0, 1, 5, 3 ,4
	};
	float solidAngles[NUMBEROFVOXELSIDES];
	float3 surfacePoints[NUMBEROFVOXELSIDES];
	float totalSolidAngle = 0;
	for (uint i = 0; i < NUMBEROFVOXELSIDES; i++)
	{
		float3 N = normals[i];

		float3 T = normals[t_indices[i]];
		float3 B = normals[b_indices[i]];

		float r0 = RandomFloat(seedptr);
		float r1 = RandomFloat(seedptr);
		float3 surfacePoint = center + (0.5 * size * N) - (0.5 * size * T) - (0.5 * size * B) + (r0 * size * T) + (r1 * size * B);

		float3 R = surfacePoint - shadingPoint;
		float3 L = normalize(R);
		float distsquared = max(10e-6, dot(R, R));
		float NlightDotL = dot(N, -L);
		float area = size * size;
		float solidAngle = NlightDotL * area / distsquared;

		surfacePoints[i] = surfacePoint;
		solidAngles[i] = max(0.0, solidAngle);
		totalSolidAngle += max(0.0, solidAngle);
	}

	float pdf[NUMBEROFVOXELSIDES];
	float cdf[NUMBEROFVOXELSIDES];
	float r2 = RandomFloat(seedptr);
	uint index = 0;
	pdf[0] = solidAngles[0] / totalSolidAngle;
	cdf[0] = pdf[0];
	if (cdf[0] < r2)
	{
		for (uint i = 1; i < NUMBEROFVOXELSIDES; i++)
		{
			pdf[i] = solidAngles[i] / totalSolidAngle;
			cdf[i] = cdf[i - 1] + pdf[i];
			index = i;
			if (cdf[i] >= r2)
			{
				break;
			}
		}
	}

	// pdf of the angle from the shading point to the emitting surface
	float weight = 1.0 / pdf[index];
	if (isinf(weight) > 0.0 || isnan(weight) > 0.0)
	{
		// we failed to calculate correct solid angles somehow and this rarily happens but does happen.
		weight = 0.0;
	}
	// the pdf of the sampled square emitting surfaces
	*invProbability = weight;

	*Nlight = normals[index];
	return surfacePoints[index];
}

float3 lightContribution(const uint voxelvalue, const uint size)
{
	float3 value = ToFloatRGB(voxelvalue) * EmitStrength(voxelvalue);
#if VOXELSAREPOINTLIGHTS
	// https://www.pbr-book.org/3ed-2018/Light_Sources/Point_Lights#fragment-PointLightMethodDefinitions-1
	value *= 4 * PI;
#else
	// https://www.pbr-book.org/3ed-2018/Light_Sources/Area_Lights#DiffuseAreaLight::Lemit
	float area = 6 * size * size;
	value *= area * PI;
#endif
	return value;
}

void GenerateLightCandidate(const uint numberOfLights, float* probability, uint* index, uint* seedptr)
{
	*index = getRandomIndex(numberOfLights, seedptr);
	// pdf of the light when uniformly sample is 1.0 / number of lights.
	*probability = 1.0 / numberOfLights;
}

void PointOnVoxelLight(const struct Light* light,
	float3 shadingPoint, uint* seedptr,
	float3* positionOnVoxel, float* invPositionProbability, float3* Nlight)
{
	const int sizeOfLight = light->size;
	// take middle of voxel
	uint3 emitterVoxelCoords = indexToCoordinates(light->position);
	const float3 center = convert_float3(emitterVoxelCoords) + 0.5 * (float3)(sizeOfLight);

#if VOXELSAREPOINTLIGHTS == 1
	//int innerSizeOfLight = sizeOfLight - 2;
	//int numberOfOutsideVoxels = sizeOfLight * sizeOfLight * sizeOfLight - innerSizeOfLight * innerSizeOfLight * innerSizeOfLight;
	// the lights are all cube shaped so voxels that are not on the outside will never contribute and therefore have pdf = 0
	//*invPositionProbability = max(0.0, convert_float(numberOfOutsideVoxels));
	* invPositionProbability = 1;
	*positionOnVoxel = center;
	*Nlight = (float3)(0);
#else
	* positionOnVoxel = RandomPointOnVoxel(center, shadingPoint, sizeOfLight,
		seedptr, invPositionProbability, Nlight);
#endif
}
