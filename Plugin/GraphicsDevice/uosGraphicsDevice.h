#include "pch.h"

#if defined(uosSupportTextureMesh)

// Graphics device identifiers in Unity
enum GfxDeviceRenderer
{
    kGfxRendererOpenGL = 0,          // OpenGL
    kGfxRendererD3D9,                // Direct3D 9
    kGfxRendererD3D11,               // Direct3D 11
    kGfxRendererGCM,                 // Sony PlayStation 3 GCM
    kGfxRendererNull,                // "null" device (used in batch mode)
    kGfxRendererHollywood,           // Nintendo Wii
    kGfxRendererXenon,               // Xbox 360
    kGfxRendererOpenGLES,            // OpenGL ES 1.1
    kGfxRendererOpenGLES20Mobile,    // OpenGL ES 2.0 mobile variant
    kGfxRendererMolehill,            // Flash 11 Stage3D
    kGfxRendererOpenGLES20Desktop,   // OpenGL ES 2.0 desktop variant (i.e. NaCl)
    kGfxRendererCount
};

// Event types for UnitySetGraphicsDevice
enum GfxDeviceEventType {
    kGfxDeviceEventInitialize = 0,
    kGfxDeviceEventShutdown,
    kGfxDeviceEventBeforeReset,
    kGfxDeviceEventAfterReset,
};

enum uosETextureFormat
{
    uosE_ARGB32 = 0,
    uosE_Depth = 1,
    uosE_ARGBHalf = 2,
    uosE_Shadowmap = 3,
    uosE_RGB565 = 4,
    uosE_ARGB4444 = 5,
    uosE_ARGB1555 = 6,
    uosE_Default = 7,
    uosE_ARGB2101010 = 8,
    uosE_DefaultHDR = 9,
    uosE_ARGBFloat = 11,
    uosE_RGFloat = 12,
    uosE_RGHalf = 13,
    uosE_RFloat = 14,
    uosE_RHalf = 15,
    uosE_R8 = 16,
    uosE_ARGBInt = 17,
    uosE_RGInt = 18,
    uosE_RInt = 19,
};

class uosIGraphicsDevice
{
public:
    virtual ~uosIGraphicsDevice() {}
    virtual void* getDevicePtr() = 0;
    virtual int getDeviceType() = 0;
    virtual bool readTexture(void *o_buf, size_t bufsize, void *tex, int width, int height, uosETextureFormat format) = 0;
    virtual bool writeTexture(void *o_tex, int width, int height, uosETextureFormat format, const void *buf, size_t bufsize) = 0;
};
uosCLinkage uosExport uosIGraphicsDevice* aiGetGraphicsDevice();
int uosGetPixelSize(uosETextureFormat format);


template<class IntType>
inline IntType ceildiv(IntType a, IntType b)
{
    return a / b + (a%b == 0 ? 0 : 1);
}

#endif
