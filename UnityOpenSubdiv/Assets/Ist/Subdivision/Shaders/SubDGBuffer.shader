Shader "Subdivision/SubDGBuffer" {
Properties {
    _Color ("Color", Color) = (1,1,1,1)
    _MainTex ("Albedo (RGB)", 2D) = "white" {}
    _Glossiness ("Smoothness", Range(0,1)) = 0.5
    _Metallic ("Metallic", Range(0,1)) = 0.0
    _EmissionMap ("Emission Map", 2D) = "white" {}
    _Emission ("EmissionColor", Color) = (1,1,1,1)
}

CGINCLUDE
#define ENABLE_INSTANCE_BUFFER 1

#include "UnityStandardCore.cginc"
#include "Assets/Ist/Foundation/Shaders/Math.cginc"
#include "Assets/Ist/Foundation/Shaders/Geometry.cginc"
#include "Assets/Ist/Foundation/Shaders/BezierPatch.cginc"
#include "Assets/Ist/BatchRenderer/Shaders/BatchRenderer.cginc"

sampler2D _NormalMap;
sampler2D _SpecularMap;
sampler2D _GrossMap;
//fixed4 _Color; // already defined in UnityCG
fixed4 _Emission;

StructuredBuffer<BezierPatch> _BPatches;





struct ia_out
{
    uint vertexID : SV_VertexID;
    uint instanceID : SV_InstanceID;
};

struct vs_out
{
    float4 vertex : SV_POSITION;
    float3 normal : TEXCOORD0;
    float4 tangent : TEXCOORD1;
    float2 texcoord : TEXCOORD2;
};

struct ps_out
{
    half4 diffuse           : SV_Target0; // RT0: diffuse color (rgb), occlusion (a)
    half4 spec_smoothness   : SV_Target1; // RT1: spec color (rgb), smoothness (a)
    half4 normal            : SV_Target2; // RT2: normal (rgb), --unused, very low precision-- (a) 
    half4 emission          : SV_Target3; // RT3: emission (rgb), --unused-- (a)
    float depth             : SV_Depth;
};


vs_out vert(ia_out I)
{
    int vid = I.vertexID;
    int iid = I.instanceID;

    vs_out O;
    O.vertex    = 0.0;
    O.normal    = 0.0;
    O.tangent   = 0.0;
    O.texcoord  = 0.0;
    return O;
}


ps_out frag(vs_out I)
{
    ps_out O;
    O.diffuse = tex2D(_MainTex, I.texcoord) * _Color;
    O.spec_smoothness = _Glossiness;
    O.normal = float4(I.normal*0.5+0.5, 0.0);
    O.emission = tex2D(_EmissionMap, I.texcoord) * _Emission;
    O.depth = 0.0;

    O.diffuse.xyz = _BPatches[0].Evaluate(I.texcoord);
    return O;
}
ENDCG

SubShader {
    Tags { "RenderType"="Opaque" }
    Cull Off

    Pass {
        Name "DEFERRED"
        Tags { "LightMode" = "Deferred" }
        Stencil {
            Comp Always
            Pass Replace
            Ref 128
        }
CGPROGRAM
#pragma target 5.0
#pragma vertex vert
#pragma fragment frag
ENDCG
    }
}
Fallback Off
}
