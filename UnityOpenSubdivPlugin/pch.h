#include <algorithm>
#include <map>
#include <memory>
#include <thread>
#include <mutex>
#include <functional>
#include <atomic>


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
#include <opensubdiv/far/primvarRefiner.h>

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

#ifndef aiNoAutoLink
#endif // aiNoAutoLink

#ifdef uosSupportD3D11
#include <d3d11.h>
#endif // uosSupportD3D11

#endif // uosWindows

