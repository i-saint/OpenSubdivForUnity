#include <algorithm>
#include <map>
#include <memory>
#include <thread>
#include <mutex>
#include <functional>
#include <atomic>

#include <glm/glm.hpp>
typedef glm::vec2 float2;
typedef glm::vec3 float3;
typedef glm::vec4 float4;
typedef glm::mat2x2 float2x2;
typedef glm::mat3x3 float3x3;
typedef glm::mat4x4 float4x4;


#ifdef _MSC_VER
    #define and &&
    #define and_eq &=
    #define bitand &
    #define bitor |
    #define compl ~
    #define not !
    #define not_eq !=
    #define or ||
    #define or_eq |=
    #define xor ^
    #define xor_eq ^=
#endif

#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/topologyRefiner.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/patchParam.h>
#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/stencilTableFactory.h>
#include <opensubdiv/osd/bufferDescriptor.h>
#include <opensubdiv/osd/cpuVertexBuffer.h>
#include <opensubdiv/osd/cpuEvaluator.h>
#include <opensubdiv/regression/common/far_utils.h>
#include <opensubdiv/regression/common/shape_utils.h>

#ifdef _WIN32
#define uosWindows
#endif // _WIN32

#define uosCLinkage extern "C"
#ifdef _MSC_VER
#define uosExport __declspec(dllexport)
#else
#define uosExport __attribute__((visibility("default")))
#endif

#ifdef uosDebug
void uosDebugLogImpl(const char* fmt, ...);
#define uosDebugLog(...) uosDebugLogImpl(__VA_ARGS__)
#ifdef uosVerboseDebug
#define uosDebugLogVerbose(...) uosDebugLogImpl(__VA_ARGS__)
#else
#define uosDebugLogVerbose(...)
#endif
#else
#define uosDebugLog(...)
#define uosDebugLogVerbose(...)
#endif


#ifdef uosWindows
#include <windows.h>
#ifdef min
#undef min
#undef max
#endif

#ifndef aiNoAutoLink
#endif // aiNoAutoLink

#ifdef uosSupportD3D11
#include <d3d11.h>
#endif // uosSupportD3D11

#endif // uosWindows

