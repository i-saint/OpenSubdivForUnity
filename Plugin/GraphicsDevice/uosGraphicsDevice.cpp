#include "pch.h"

#if defined(uosSupportTextureMesh)

#include "uosGraphicsDevice.h"


int uosGetPixelSize(uosETextureFormat format)
{
    switch (format)
    {
    case uosE_ARGB32:    return 4;

    case uosE_ARGBHalf:  return 8;
    case uosE_RGHalf:    return 4;
    case uosE_RHalf:     return 2;

    case uosE_ARGBFloat: return 16;
    case uosE_RGFloat:   return 8;
    case uosE_RFloat:    return 4;

    case uosE_ARGBInt:   return 16;
    case uosE_RGInt:     return 8;
    case uosE_RInt:      return 4;
    }
    return 0;
}


uosIGraphicsDevice* uosCreateGraphicsDeviceOpenGL(void *device);
uosIGraphicsDevice* uosCreateGraphicsDeviceD3D9(void *device);
uosIGraphicsDevice* uosCreateGraphicsDeviceD3D11(void *device);


uosIGraphicsDevice *g_the_graphics_device;
uosCLinkage uosExport uosIGraphicsDevice* aiGetGraphicsDevice() { return g_the_graphics_device; }
typedef uosIGraphicsDevice* (*aiGetGraphicsDeviceT)();


uosCLinkage uosExport void UnitySetGraphicsDevice(void* device, int deviceType, int eventType)
{
    if (eventType == kGfxDeviceEventInitialize) {
#ifdef aiSupportD3D9
        if (deviceType == kGfxRendererD3D9)
        {
            g_the_graphics_device = uosCreateGraphicsDeviceD3D9(device);
        }
#endif // aiSupportD3D9
#ifdef uosSupportD3D11
        if (deviceType == kGfxRendererD3D11)
        {
            g_the_graphics_device = uosCreateGraphicsDeviceD3D11(device);
        }
#endif // uosSupportD3D11
#ifdef uosSupportOpenGL
        if (deviceType == kGfxRendererOpenGL)
        {
            g_the_graphics_device = uosCreateGraphicsDeviceOpenGL(device);
        }
#endif // uosSupportOpenGL
    }

    if (eventType == kGfxDeviceEventShutdown) {
        delete g_the_graphics_device;
        g_the_graphics_device = nullptr;
    }
}

uosCLinkage uosExport void UnityRenderEvent(int eventID)
{
}


#ifdef uosSupportOpenGL
uosCLinkage uosExport void aiInitializeOpenGL()
{
    UnitySetGraphicsDevice(nullptr, kGfxRendererOpenGL, kGfxDeviceEventInitialize);
}
#endif

#ifdef aiSupportD3D9
uosCLinkage uosExport void aiInitializeD3D9(void *device)
{
    UnitySetGraphicsDevice(device, kGfxRendererD3D9, kGfxDeviceEventInitialize);
}
#endif

#ifdef uosSupportD3D11
uosCLinkage uosExport void aiInitializeD3D11(void *device)
{
    UnitySetGraphicsDevice(device, kGfxRendererD3D11, kGfxDeviceEventInitialize);
}
#endif

uosCLinkage uosExport void aiFinalizeGraphicsDevice()
{
    UnitySetGraphicsDevice(nullptr, kGfxRendererNull, kGfxDeviceEventShutdown);
}



#if !defined(aiMaster) && defined(uosWindows)

// PatchLibrary で突っ込まれたモジュールは UnitySetGraphicsDevice() が呼ばれないので、
// DLL_PROCESS_ATTACH のタイミングで先にロードされているモジュールからデバイスをもらって同等の処理を行う。
BOOL WINAPI DllMain(HINSTANCE module_handle, DWORD reason_for_call, LPVOID reserved)
{
    if (reason_for_call == DLL_PROCESS_ATTACH)
    {
        HMODULE m = ::GetModuleHandleA("AlembicImporter.dll");
        if (m) {
            auto proc = (aiGetGraphicsDeviceT)::GetProcAddress(m, "aiGetGraphicsDevice");
            if (proc) {
                uosIGraphicsDevice *dev = proc();
                if (dev) {
                    UnitySetGraphicsDevice(dev->getDevicePtr(), dev->getDeviceType(), kGfxDeviceEventInitialize);
                }
            }
        }
    }
    else if (reason_for_call == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}

// "DllMain already defined in MSVCRT.lib" 対策
#ifdef _X86_
extern "C" { int _afxForceUSRDLL; }
#else
extern "C" { int __afxForceUSRDLL; }
#endif

#endif

#endif 
