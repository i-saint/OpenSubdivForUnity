Shader "Subdivision/Bezier Patch Raytracer" {
Properties {
    _Color("Color", Color) = (0.5, 0.5, 0.5, 0.5)
    _SpecularColor("Specular Color", Color) = (0.5, 0.5, 0.5, 0.5)
    _EmissionColor("Emission Color", Color) = (0.0, 0.0, 0.0, 0.0)
    _MainTex ("Albedo (RGB)", 2D) = "white" {}
    _Glossiness("Smoothness", Range(0,1)) = 0.5
    _Metallic("Metallic", Range(0, 1)) = 0.0
}

CGINCLUDE
#include "UnityStandardCore.cginc"
#include "Assets/Ist/Foundation/Shaders/Math.cginc"
#include "Assets/Ist/Foundation/Shaders/Geometry.cginc"
#include "Assets/Ist/Foundation/Shaders/BuiltinVariablesExt.cginc"
#include "Assets/Ist/Foundation/Shaders/BezierPatch.cginc"
#include "Assets/Ist/Foundation/Shaders/BezierPatchIntersection.cginc"

float4 _SpecularColor;



struct ia_out
{
    float4 vertex : POSITION;
};

struct vs_out
{
    float4 vertex : SV_POSITION;
    float4 screen_pos : TEXspos0;
    float4 world_pos : TEXspos1;
};


vs_out vert(ia_out I)
{
    vs_out O;
    O.vertex = mul(UNITY_MATRIX_MVP, I.vertex);
    //O.screen_pos = ComputeScreenPos(O.vertex);
    O.screen_pos = O.vertex;
    O.world_pos = mul(_Object2World, I.vertex);
    return O;
}


float3 localize(float3 p)
{
    p = mul(_World2Object, float4(p, 1.0)).xyz;
    return p;
}

void raymarching(float2 pos, inout float distance)
{

}



struct gbuffer_out
{
    half4 diffuse           : SV_Target0; // RT0: diffuse color (rgb), occlusion (a)
    half4 spec_smoothness   : SV_Target1; // RT1: spec color (rgb), smoothness (a)
    half4 normal            : SV_Target2; // RT2: normal (rgb), --unused, very low precision-- (a) 
    half4 emission          : SV_Target3; // RT3: emission (rgb), --unused-- (a)
    float depth             : SV_Depth;
};


gbuffer_out frag_gbuffer(vs_out I)
{
    float3 world_pos = I.world_pos.xyz;
    I.screen_pos.xy /= I.screen_pos.w;
#if UNITY_UV_STARTS_AT_TOP
    I.screen_pos.y *= -1.0;
#endif
    float2 spos = I.screen_pos.xy;
    spos.x *= _ScreenParams.x / _ScreenParams.y;

    float max_distance = _ProjectionParams.z - _ProjectionParams.y;
    float distance = 0.0;
    float zmin = 0.0;          // todo 
    float zmax = max_distance;  // 

    BezierPatch bpatch;
    BezierPatchHit hit;
    Ray ray;
    UNITY_INITIALIZE_OUTPUT(BezierPatch, bpatch);


    {
        float3 cam_pos = GetCameraPosition();
        float3 cam_forward = GetCameraForward();
        float3 cam_up = GetCameraUp();
        float3 cam_right = GetCameraRight();
        float  cam_focal_len = GetCameraFocalLength();

        ray.direction = normalize(cam_right*spos.x + cam_up*spos.y + cam_forward*cam_focal_len);
        ray.origin = cam_pos + ray.direction * distance;
    }

    if (!BPIRaycast(bpatch, ray, zmin, zmax, hit)) {
        discard;
    }

    float3 pos = ray.origin + ray.direction * hit.t;
    float3 normal = BPEvaluateNormal(bpatch, float2(hit.u, hit.v));


    gbuffer_out O;
    O.diffuse = _Color;
    O.spec_smoothness = float4(_SpecularColor.rgb, _Glossiness);
    O.normal = float4(normal*0.5+0.5, 1.0);
    O.emission = _EmissionColor;
#ifndef UNITY_HDR_ON
    O.emission = exp2(-O.emission);
#endif
    O.depth = ComputeDepth(mul(UNITY_MATRIX_VP, float4(pos, 1.0)));
    return O;
}

ENDCG

SubShader {
    Tags{ "RenderType" = "Opaque" "DisableBatching" = "True" "Queue" = "Geometry+10" }
    Cull Off

    Pass {
        Tags { "LightMode" = "Deferred" }
        Stencil {
            Comp Always
            Pass Replace
            Ref 128
        }
CGPROGRAM
#pragma target 3.0
#pragma vertex vert
#pragma fragment frag_gbuffer
#pragma multi_compile ___ UNITY_HDR_ON
ENDCG
    }
}
Fallback Off
}
