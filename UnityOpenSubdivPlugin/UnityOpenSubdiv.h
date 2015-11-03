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



uosCLinkage uosExport cfloat3 BezierPatchEvaluate(const BezierPatch *bp, const float2 *uv);
uosCLinkage uosExport cfloat3 BezierPatchEvaluateNormal(const BezierPatch *bp, const float2 *uv);
uosCLinkage uosExport bool BezierPatchRayIntersection(const BezierPatch *bp, const Ray *ray, BezierPatchHit *hit);

#endif // UnityOpenSubdiv_h
