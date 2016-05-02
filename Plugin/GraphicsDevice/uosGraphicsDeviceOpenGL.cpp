#include "pch.h"

#if defined(uosSupportTextureMesh) && defined(uosSupportOpenGL)
#include "uosGraphicsDevice.h"

#define GLEW_STATIC
#include <GL/glew.h>
#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glew32.lib")


class uosGraphicsDeviceOpenGL : public uosIGraphicsDevice
{
public:
    uosGraphicsDeviceOpenGL(void *device);
    ~uosGraphicsDeviceOpenGL();
    void* getDevicePtr() override;
    int getDeviceType() override;
    bool readTexture(void *o_buf, size_t bufsize, void *tex, int width, int height, uosETextureFormat format) override;
    bool writeTexture(void *o_tex, int width, int height, uosETextureFormat format, const void *buf, size_t bufsize) override;

private:
    void *m_device;
};


uosIGraphicsDevice* uosCreateGraphicsDeviceOpenGL(void *device)
{
    return new uosGraphicsDeviceOpenGL(device);
}


void* uosGraphicsDeviceOpenGL::getDevicePtr() { return m_device; }
int uosGraphicsDeviceOpenGL::getDeviceType() { return kGfxRendererOpenGL; }

uosGraphicsDeviceOpenGL::uosGraphicsDeviceOpenGL(void *device)
    : m_device(device)
{
    glewInit();
}

uosGraphicsDeviceOpenGL::~uosGraphicsDeviceOpenGL()
{
}


static void uosGetInternalFormatOpenGL(uosETextureFormat format, GLenum &o_fmt, GLenum &o_type)
{
    switch (format)
    {
    case uosE_ARGB32:    o_fmt = GL_RGBA; o_type = GL_UNSIGNED_BYTE; return;

    case uosE_ARGBHalf:  o_fmt = GL_RGBA; o_type = GL_HALF_FLOAT; return;
    case uosE_RGHalf:    o_fmt = GL_RG; o_type = GL_HALF_FLOAT; return;
    case uosE_RHalf:     o_fmt = GL_RED; o_type = GL_HALF_FLOAT; return;

    case uosE_ARGBFloat: o_fmt = GL_RGBA; o_type = GL_FLOAT; return;
    case uosE_RGFloat:   o_fmt = GL_RG; o_type = GL_FLOAT; return;
    case uosE_RFloat:    o_fmt = GL_RED; o_type = GL_FLOAT; return;

    case uosE_ARGBInt:   o_fmt = GL_RGBA_INTEGER; o_type = GL_INT; return;
    case uosE_RGInt:     o_fmt = GL_RG_INTEGER; o_type = GL_INT; return;
    case uosE_RInt:      o_fmt = GL_RED_INTEGER; o_type = GL_INT; return;
    }
}

bool uosGraphicsDeviceOpenGL::readTexture(void *o_buf, size_t bufsize, void *tex, int width, int height, uosETextureFormat format)
{
    GLenum internal_format = 0;
    GLenum internal_type = 0;
    uosGetInternalFormatOpenGL(format, internal_format, internal_type);

    //// glGetTextureImage() is available only OpenGL 4.5 or later...
    // glGetTextureImage((GLuint)(size_t)tex, 0, internal_format, internal_type, bufsize, o_buf);

    glFinish();
    glBindTexture(GL_TEXTURE_2D, (GLuint)(size_t)tex);
    glGetTexImage(GL_TEXTURE_2D, 0, internal_format, internal_type, o_buf);
    glBindTexture(GL_TEXTURE_2D, 0);
    return true;
}

bool uosGraphicsDeviceOpenGL::writeTexture(void *o_tex, int width, int height, uosETextureFormat format, const void *buf, size_t bufsize)
{
    GLenum internal_format = 0;
    GLenum internal_type = 0;
    uosGetInternalFormatOpenGL(format, internal_format, internal_type);

    //// glTextureSubImage2D() is available only OpenGL 4.5 or later...
    // glTextureSubImage2D((GLuint)(size_t)o_tex, 0, 0, 0, width, height, internal_format, internal_type, buf);

    glBindTexture(GL_TEXTURE_2D, (GLuint)(size_t)o_tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, internal_format, internal_type, buf);
    glBindTexture(GL_TEXTURE_2D, 0);
    return true;
}

#endif // uosSupportOpenGL
