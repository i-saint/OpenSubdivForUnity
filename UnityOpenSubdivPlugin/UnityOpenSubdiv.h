#ifndef UnityOpenSubdiv_h
#define UnityOpenSubdiv_h

// options:
// graphics device options are relevant only if uosSupportTextureMesh is defined
// 
// #define uosSupportTextureMesh
//  #define uosSupportD3D11
//  #define uosSupportOpenGL
// 
// #define uosWithTBB


#include "pch.h"


class uosContext;

uosCLinkage uosExport uosContext*   uosCreateContext();
uosCLinkage uosExport void          uosDestroyContext(uosContext* ctx);

#endif // UnityOpenSubdiv_h
