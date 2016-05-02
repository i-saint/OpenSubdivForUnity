#include "pch.h"
#include "osuInternal.h"
#include "OpenSubdivForUnity.h"
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






osuCLinkage osuExport osuContext* osuCreateContext()
{
    osuContext *ctx = nullptr;
    osuDebugLog("osuCreateContext(): %p\n", ctx);
    return ctx;
}

osuCLinkage osuExport void osuDestroyContext(osuContext* ctx)
{
    osuDebugLog("osuDestroyContext(): %p\n", ctx);
}


osuCLinkage osuExport void osuBPGetMinMax(const BezierPatch *bp, float3 *omin, float3 *omax, float eps)
{
    bp->GetMinMax(*omin, *omax, eps);
}

osuCLinkage osuExport cfloat3 osuBPGetRoughNormal(const BezierPatch *bp)
{
    return (cfloat3&)bp->GetRoughNormal();
}

osuCLinkage osuExport cfloat3 osuBPEvaluate(const BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->Evaluate(*uv);
}
osuCLinkage osuExport cfloat3 osuBPEvaluateNormal(const BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->EvaluateNormal(*uv);
}

osuCLinkage osuExport void osuBPSplit(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, BezierPatch *dst2, BezierPatch *dst3, const float2 *uv)
{
    bp->Split(*dst0, *dst1, *dst2, *dst3, *uv);
}
osuCLinkage osuExport void osuBPSplitU(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float u)
{
    bp->SplitU(*dst0, *dst1, u);
}
osuCLinkage osuExport void osuBPSplitV(const BezierPatch *bp, BezierPatch *dst0, BezierPatch *dst1, float v)
{
    bp->SplitV(*dst0, *dst1, v);
}
osuCLinkage osuExport void osuBPCrop(const BezierPatch *bp, BezierPatch *dst, const float2 *uv0, const float2 *uv1)
{
    bp->Crop(*dst, *uv0, *uv1);
}
osuCLinkage osuExport void osuBPCropU(const BezierPatch *bp, BezierPatch *dst, float u0, float u1)
{
    bp->CropU(*dst, u0, u1);
}
osuCLinkage osuExport void osuBPCropV(const BezierPatch *bp, BezierPatch *dst, float v0, float v1)
{
    bp->CropV(*dst, v0, v1);
}

osuCLinkage osuExport void osuBPTransform(BezierPatch *bp, float4x4 *mat)
{
    bp->Transform(*mat);
}

osuCLinkage osuExport bool osuBPRaycast(const BezierPatch *bp, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit)
{
    Ray ray = { *orig, *dir };
    return BezierPatchRaycast(*bp, ray, zmin, zmax, epsilon, *hit);
}

osuCLinkage osuExport bool osuBPRaycastWithTransform(const BezierPatch *bp, const float4x4 *trans, const float3 *orig, const float3 *dir, float zmin, float zmax, float epsilon, BezierPatchHit *hit)
{
    Ray ray = { *orig, *dir };
    return BezierPatchRaycast(*bp, *trans, ray, zmin, zmax, epsilon, *hit);
}
