#pragma once

#include "BezierPatch.h"

namespace osdi {

struct Ray
{
    float3 orig;
    float3 dir;
};

struct BezierPatchHit
{
    float2 uv;
    float t;
    uint32_t clip_level;
};

bool BezierPatchRaycast(const BezierPatch &bp, const Ray &ray, float zmin, float zmax, float epsilon, BezierPatchHit &hit);
bool BezierPatchRaycast(const BezierPatch &bp, const float4x4 &trans, const Ray &ray, float zmin, float zmax, float epsilon, BezierPatchHit &hit);

} // namespace osdi
