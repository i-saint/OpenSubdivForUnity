#include "pch.h"
#include "UnityOpenSubdiv.h"
#include "BezierPatch/BezierPatch.h"
#include "BezierPatch/BezierPatchIntersection.h"

#ifdef uosWindows
    #include <windows.h>
    //#pragma comment(lib, "osdCPU.lib")
    //#pragma comment(lib, "osdGPU.lib")
#else // uosWindows
#endif // uosWindows


#ifdef uosDebug
void uosDebugLogImpl(const char* fmt, ...)
{
    va_list vl;
    va_start(vl, fmt);

#ifdef uosWindows
    char buf[2048];
    vsprintf(buf, fmt, vl);
    ::OutputDebugStringA(buf);
#else // uosWindows
    vprintf(fmt, vl);
#endif // uosWindows

    va_end(vl);
}
#endif // uosDebug






uosCLinkage uosExport uosContext* uosCreateContext()
{
    uosContext *ctx = nullptr;
    uosDebugLog("uosCreateContext(): %p\n", ctx);
    return ctx;
}

uosCLinkage uosExport void uosDestroyContext(uosContext* ctx)
{
    uosDebugLog("uosDestroyContext(): %p\n", ctx);
}


uosCLinkage uosExport void uosBPGetMinMax(const BezierPatch *bp, float3 *omin, float3 *omax, float eps)
{
    bp->GetMinMax(*omin, *omax, eps);
}

uosCLinkage uosExport cfloat3 uosBPGetRoughNormal(const BezierPatch *bp)
{
    return (cfloat3&)bp->GetRoughNormal();
}

uosCLinkage uosExport cfloat3 uosBPEvaluate(const BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->Evaluate(*uv);
}
uosCLinkage uosExport cfloat3 uosBPEvaluateNormal(const BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->EvaluateNormal(*uv);
}

uosCLinkage uosExport void uosBPSplit(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, BezierPatch *dst2, BezierPatch *dst3, const float2 *uv)
{
    bp->Split(*dst0, *dst1, *dst2, *dst3, *uv);
}
uosCLinkage uosExport void uosBPSplitU(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float u)
{
    bp->SplitU(*dst0, *dst1, u);
}
uosCLinkage uosExport void uosBPSplitV(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float v)
{
    bp->SplitV(*dst0, *dst1, v);
}
uosCLinkage uosExport void uosBPCrop(const BezierPatch *bp, BezierPatch *dst, const float2 *uv0, const float2 *uv1)
{
    bp->Crop(*dst, *uv0, *uv1);
}
uosCLinkage uosExport void uosBPCropU(const BezierPatch *bp, BezierPatch *dst, float u0, float u1)
{
    bp->CropU(*dst, u0, u1);
}
uosCLinkage uosExport void uosBPCropV(const BezierPatch *bp, BezierPatch *dst, float v0, float v1)
{
    bp->CropV(*dst, v0, v1);
}

uosCLinkage uosExport void uosBPTransform(BezierPatch *bp, float4x4 *mat)
{
    bp->Transform(*mat);
}

uosCLinkage uosExport bool uosBPRaycast(const BezierPatch *bp, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit)
{
    Ray ray = { *orig, *dir };
    return BezierPatchRaycast(*bp, ray, zmin, zmax, epsilon, *hit);
}

uosCLinkage uosExport bool uosBPRaycastWithTransform(const BezierPatch *bp, const float4x4 *trans, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit)
{
    Ray ray = { *orig, *dir };
    return BezierPatchRaycast(*bp, *trans, ray, zmin, zmax, epsilon, *hit);
}
