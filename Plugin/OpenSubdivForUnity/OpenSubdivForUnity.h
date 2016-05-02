#pragma once

#define osuCLinkage extern "C"
#ifdef _MSC_VER
    #define osuExport __declspec(dllexport)
#else
    #define osuExport __attribute__((visibility("default")))
#endif

class osuContext;
struct BezierPatch;
struct BezierPatchHit;
struct Ray;

struct cfloat2 { float v[2]; }; // to shut up compiler
struct cfloat3 { float v[3]; }; //
struct cfloat4 { float v[4]; }; //

osuCLinkage osuExport osuContext*   osuCreateContext();
osuCLinkage osuExport void          osuDestroyContext(osuContext* ctx);



osuCLinkage osuExport void osuBPGetMinMax(const BezierPatch *bp, float3 *omin, float3 *omax, float eps);
osuCLinkage osuExport cfloat3 osuBPGetRoughNormal(const BezierPatch *bp);
osuCLinkage osuExport cfloat3 osuBPEvaluate(const BezierPatch *bp, const float2 *uv);
osuCLinkage osuExport cfloat3 osuBPEvaluateNormal(const BezierPatch *bp, const float2 *uv);
osuCLinkage osuExport void osuBPSplit(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, BezierPatch *dst2, BezierPatch *dst3, const float2 *uv);
osuCLinkage osuExport void osuBPSplitU(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float u);
osuCLinkage osuExport void osuBPSplitV(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float v);
osuCLinkage osuExport void osuBPCrop(const BezierPatch *bp, BezierPatch *dst, const float2 *uv0, const float2 *uv1);
osuCLinkage osuExport void osuBPCropU(const BezierPatch *bp, BezierPatch *dst, float u0, float u1);
osuCLinkage osuExport void osuBPCropV(const BezierPatch *bp, BezierPatch *dst, float v0, float v1);
osuCLinkage osuExport void osuBPTransform(BezierPatch *bp, float4x4 *mat);

osuCLinkage osuExport bool osuBPRaycast(const BezierPatch *bp, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit);
osuCLinkage osuExport bool osuBPRaycastWithTransform(const BezierPatch *bp, const float4x4 *trans, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit);
