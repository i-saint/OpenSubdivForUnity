#pragma once

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

typedef glm::vec2 float2;
typedef glm::vec3 float3;
typedef glm::vec4 float4;
typedef glm::mat2x2 float2x2;
typedef glm::mat3x3 float3x3;
typedef glm::mat4x4 float4x4;
