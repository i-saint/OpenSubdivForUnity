#ifndef IstBezierPatch_h
#define IstBezierPatch_h

struct BezierPatch
{
    float4 m_cp[16];


    float3 Evaluate(float2 uv);
    float3 EvaluateDu(float2 uv);
    float3 EvaluateDv(float2 uv);
    float3 EvaluateNormal(float2 uv);

    // impl
    float3 EvaluateImpl(float t, float3 p[4]);
    float3 EvaluateDImpl(float t, float3 p[4]);
};



float3 BezierPatch::Evaluate(float2 uv)
{
    float3 b[4];
    for (int i = 0; i < 4; ++i) {
        float3 cp[4] = {
            m_cp[i * 4 + 0].xyz,
            m_cp[i * 4 + 1].xyz,
            m_cp[i * 4 + 2].xyz,
            m_cp[i * 4 + 3].xyz
        };
        b[i] = EvaluateImpl(uv.x, cp);
    }
    return EvaluateImpl(uv.y, b);
}

float3 BezierPatch::EvaluateDu(float2 uv)
{
    float3 b[4];
    for (int i = 0; i < 4; ++i) {
        float3 cp[4] = {
            m_cp[i * 4 + 0].xyz,
            m_cp[i * 4 + 1].xyz,
            m_cp[i * 4 + 2].xyz,
            m_cp[i * 4 + 3].xyz
        };
        b[i] = EvaluateDImpl(uv.x, cp);
    }
    return EvaluateImpl(uv.y, b);
}

float3 BezierPatch::EvaluateDv(float2 uv)
{
    float3 b[4];
    for (int i = 0; i < 4; ++i) {
        float3 cp[4] = {
            m_cp[i * 4 + 0].xyz,
            m_cp[i * 4 + 1].xyz,
            m_cp[i * 4 + 2].xyz,
            m_cp[i * 4 + 3].xyz
        };
        b[i] = EvaluateImpl(uv.x, cp);
    }
    return EvaluateDImpl(uv.y, b);
}

float3 BezierPatch::EvaluateNormal(float2 uv)
{
    float3 du = EvaluateDu(uv);
    float3 dv = EvaluateDv(uv);
    return normalize(cross(du, dv));
}


float3 BezierPatch::EvaluateImpl(float t, float3 cp[4])
{
    float it = 1.0f - t;
    return cp[0] * (       (it*it*it)          )
         + cp[1] * (3.0f * (it*it   ) * (t)    )
         + cp[2] * (3.0f * (it      ) * (t*t)  )
         + cp[3] * (                    (t*t*t));
}

float3 BezierPatch::EvaluateDImpl(float t, float3 cp[4])
{
    float t2 = t * t;
    return cp[0] * (3.0f * t2 *-1.0f + 2.0f * t * 3.0f - 3.0f)
         + cp[1] * (3.0f * t2 * 3.0f + 2.0f * t *-6.0f + 3.0f)
         + cp[2] * (3.0f * t2 *-3.0f + 2.0f * t * 3.0f       )
         + cp[3] * (3.0f * t2 * 1.0f                         );
}


#endif // IstBezierPatch_h
