#pragma once

#ifdef _MSC_VER
    #define osdiAPI extern "C" __declspec(dllexport)
#else
    #define osdiAPI extern "C" __attribute__((visibility("default")))
#endif

namespace osdi {
    class Context;
    struct BezierPatch;
    struct BezierPatchHit;
    struct Ray;
} // namespace osdi

struct cfloat2 { float v[2]; }; // to shut up compiler
struct cfloat3 { float v[3]; }; //
struct cfloat4 { float v[4]; }; //

osdiAPI osdi::Context*  osdiCreateContext();
osdiAPI void            osdiDestroyContext(osdi::Context* ctx);



osdiAPI void    osdiBPGetMinMax(const osdi::BezierPatch *bp, float3 *omin, float3 *omax, float eps);
osdiAPI cfloat3 osdiBPGetRoughNormal(const osdi::BezierPatch *bp);
osdiAPI cfloat3 osdiBPEvaluate(const osdi::BezierPatch *bp, const float2 *uv);
osdiAPI cfloat3 osdiBPEvaluateNormal(const osdi::BezierPatch *bp, const float2 *uv);
osdiAPI void    osdiBPSplit(const osdi::BezierPatch *bp, osdi::BezierPatch *dst0, osdi::BezierPatch *dst1, osdi::BezierPatch *dst2, osdi::BezierPatch *dst3, const float2 *uv);
osdiAPI void    osdiBPSplitU(const osdi::BezierPatch *bp, osdi::BezierPatch *dst0, osdi::BezierPatch *dst1, float u);
osdiAPI void    osdiBPSplitV(const osdi::BezierPatch *bp, osdi::BezierPatch *dst0, osdi::BezierPatch *dst1, float v);
osdiAPI void    osdiBPCrop(const osdi::BezierPatch *bp, osdi::BezierPatch *dst, const float2 *uv0, const float2 *uv1);
osdiAPI void    osdiBPCropU(const osdi::BezierPatch *bp, osdi::BezierPatch *dst, float u0, float u1);
osdiAPI void    osdiBPCropV(const osdi::BezierPatch *bp, osdi::BezierPatch *dst, float v0, float v1);
osdiAPI void    osdiBPTransform(osdi::BezierPatch *bp, float4x4 *mat);

osdiAPI bool    osdiBPRaycast(const osdi::BezierPatch *bp, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, osdi::BezierPatchHit *hit);
osdiAPI bool    osdiBPRaycastWithTransform(const osdi::BezierPatch *bp, const float4x4 *trans, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, osdi::BezierPatchHit *hit);
