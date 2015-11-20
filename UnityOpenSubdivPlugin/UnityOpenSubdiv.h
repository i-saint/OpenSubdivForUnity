#ifndef UnityOpenSubdiv_h
#define UnityOpenSubdiv_h

// options:
// graphics device options are relevant only if uosSupportTextureMesh is defined
// 
// #define uosSupportTextureMesh
//  #define uosSupportD3D11
//  #define uosSupportOpenGL
// 
// #define uosWithTBB


#include "pch.h"


class uosContext;
struct BezierPatch;
struct BezierPatchHit;
struct Ray;

struct cfloat2 { float v[2]; }; // to shut up compiler
struct cfloat3 { float v[3]; }; //
struct cfloat4 { float v[4]; }; //

uosCLinkage uosExport uosContext*   uosCreateContext();
uosCLinkage uosExport void          uosDestroyContext(uosContext* ctx);



uosCLinkage uosExport cfloat3 uosBPEvaluate(const BezierPatch *bp, const float2 *uv);
uosCLinkage uosExport cfloat3 uosBPEvaluateNormal(const BezierPatch *bp, const float2 *uv);

uosCLinkage uosExport void uosBPSplit(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, BezierPatch *dst2, BezierPatch *dst3, const float2 *uv);
uosCLinkage uosExport void uosBPSplitU(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float u);
uosCLinkage uosExport void uosBPSplitV(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float v);
uosCLinkage uosExport void uosBPCrop(const BezierPatch *bp, BezierPatch *dst, const float2 *uv0, const float2 *uv1);
uosCLinkage uosExport void uosBPCropU(const BezierPatch *bp, BezierPatch *dst, float u0, float u1);
uosCLinkage uosExport void uosBPCropV(const BezierPatch *bp, BezierPatch *dst, float v0, float v1);
uosCLinkage uosExport void uosBPTransform(BezierPatch *bp, float4x4 *mat);

uosCLinkage uosExport bool uosBPRaycast(const BezierPatch *bp, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit);
uosCLinkage uosExport bool uosBPRaycastWithTransform(const BezierPatch *bp, const float4x4 *trans, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit);

#endif // UnityOpenSubdiv_h
