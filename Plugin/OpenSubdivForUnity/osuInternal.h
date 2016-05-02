#pragma once
#include "BezierPatch/BezierPatch.h"
#include "BezierPatch/BezierPatchIntersection.h"

#ifdef _WIN32
    #define osuWindows
#endif // _WIN32


#ifdef osuDebug
    void osuDebugLogImpl(const char* fmt, ...);
    #define osuDebugLog(...) osuDebugLogImpl(__VA_ARGS__)
    #ifdef osuVerboseDebug
        #define osuDebugLogVerbose(...) osuDebugLogImpl(__VA_ARGS__)
    #else
        #define osuDebugLogVerbose(...)
    #endif
#else
    #define osuDebugLog(...)
    #define osuDebugLogVerbose(...)
#endif
