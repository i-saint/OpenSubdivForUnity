#include "pch.h"
#include "UnityOpenSubdiv.h"
#include "BezierPatch//BezierPatch.h"
#include "BezierPatch//Intersection.h"

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



uosCLinkage uosExport cfloat3 BezierPatchEvaluate(const BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->Evaluate(*uv);
}

uosCLinkage uosExport cfloat3 BezierPatchEvaluateNormal(const BezierPatch *bp, const float2 *uv)
{
    return (cfloat3&)bp->EvaluateNormal(*uv);
}
