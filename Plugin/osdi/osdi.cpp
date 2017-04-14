#include "pch.h"
#include "osuInternal.h"
#include "osdi.h"
#include "BezierPatch/BezierPatch.h"
#include "BezierPatch/BezierPatchIntersection.h"

#ifdef osuWindows
    #include <windows.h>
    #pragma comment(lib, "osdCPU.lib")
    #pragma comment(lib, "regression_common_obj.lib")
#endif // osuWindows


#ifdef osuDebug
void osuDebugLogImpl(const char* fmt, ...)
{
    va_list vl;
    va_start(vl, fmt);

#ifdef osuWindows
    char buf[2048];
    vsprintf(buf, fmt, vl);
    ::OutputDebugStringA(buf);
#else // osuWindows
    vprintf(fmt, vl);
#endif // osuWindows

    va_end(vl);
}
#endif // osuDebug






osdiAPI osdi::Context* osdiCreateContext()
{
    osdi::Context *ctx = nullptr;
    osuDebugLog("osuCreateContext(): %p\n", ctx);
    return ctx;
}

osdiAPI void osdiDestroyContext(osdi::Context* ctx)
{
    osuDebugLog("osuDestroyContext(): %p\n", ctx);
}


osdiAPI void osdiBPGetMinMax(const osdi::BezierPatch *bp, float3 *omin, float3 *omax, float eps)
{
    bp->GetMinMax(*omin, *omax, eps);
}

osdiAPI cfloat3 osdiBPGetRoughNormal(const osdi::BezierPatch *bp)
{
    return (cfloat3&)bp->GetRoughNormal();
}

osdiAPI cfloat3 osdiBPEvaluate(const osdi::BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->Evaluate(*uv);
}
osdiAPI cfloat3 osdiBPEvaluateNormal(const osdi::BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->EvaluateNormal(*uv);
}

osdiAPI void osdiBPSplit(const osdi::BezierPatch *bp, osdi::BezierPatch *dst0, osdi::BezierPatch *dst1, osdi::BezierPatch *dst2, osdi::BezierPatch *dst3, const float2 *uv)
{
    bp->Split(*dst0, *dst1, *dst2, *dst3, *uv);
}
osdiAPI void osdiBPSplitU(const osdi::BezierPatch *bp, osdi::BezierPatch *dst0, osdi::BezierPatch *dst1, float u)
{
    bp->SplitU(*dst0, *dst1, u);
}
osdiAPI void osdiBPSplitV(const osdi::BezierPatch *bp, osdi::BezierPatch *dst0, osdi::BezierPatch *dst1, float v)
{
    bp->SplitV(*dst0, *dst1, v);
}
osdiAPI void osdiBPCrop(const osdi::BezierPatch *bp, osdi::BezierPatch *dst, const float2 *uv0, const float2 *uv1)
{
    bp->Crop(*dst, *uv0, *uv1);
}
osdiAPI void osdiBPCropU(const osdi::BezierPatch *bp, osdi::BezierPatch *dst, float u0, float u1)
{
    bp->CropU(*dst, u0, u1);
}
osdiAPI void osdiBPCropV(const osdi::BezierPatch *bp, osdi::BezierPatch *dst, float v0, float v1)
{
    bp->CropV(*dst, v0, v1);
}

osdiAPI void osdiBPTransform(osdi::BezierPatch *bp, float4x4 *mat)
{
    bp->Transform(*mat);
}

osdiAPI bool osdiBPRaycast(const osdi::BezierPatch *bp, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, osdi::BezierPatchHit *hit)
{
    osdi::Ray ray = { *orig, *dir };
    return BezierPatchRaycast(*bp, ray, zmin, zmax, epsilon, *hit);
}

osdiAPI bool osdiBPRaycastWithTransform(const osdi::BezierPatch *bp, const float4x4 *trans, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, osdi::BezierPatchHit *hit)
{
    osdi::Ray ray = { *orig, *dir };
    return BezierPatchRaycast(*bp, *trans, ray, zmin, zmax, epsilon, *hit);
}
