#ifndef IstBezierPatch_h
#define IstBezierPatch_h

struct BezierPatch
{
    float4 m_cp[16];


    float3 Evaluate(float2 uv);
    float3 EvaluateDu(float2 uv);
    float3 EvaluateDv(float2 uv);
    float3 EvaluateNormal(float2 uv);

    void Split(BezierPatch patches[4], float u, float v);
    void SplitU(BezierPatch patches[2], float u);
    void SplitV(BezierPatch patches[2], float v);
    void CropU(BezierPatch patch, float u0, float u1);
    void CropV(BezierPatch patch, float v0, float v1);
    bool Crop(BezierPatch patch, float u0, float u1, float v0, float v1);
    float3 GetLv();
    float3 GetLu();

    void Transpose();
    void Transform(float3x3 m);
    void Transform(float4x4 m);

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

/*
void BezierPatch::Split(BezierPatch patches[4], float u, float v)
{
    BezierPatch tmp0, tmp1;

    // split U
    for (int i = 0; i < N; ++i) {
        bezierSplit(&tmp0[i*N + 0], &tmp1[i*N + 0], &m_cp[i*N + 0], u);
    }

    // uv -> vu
    tmp0.Transpose(); // 00 01
    tmp1.Transpose(); // 10 11

                      // split V
    for (int i = 0; i < N; ++i) {
        bezierSplit(&patches[0].m_cp[i*N + 0],
            &patches[2].m_cp[i*N + 0],
            &tmp0[i*N + 0], v);
        bezierSplit(&patches[1].m_cp[i*N + 0],
            &patches[3].m_cp[i*N + 0],
            &tmp1[i*N + 0], v);
    }

    // vu -> uv
    patches[0].Transpose(); //00
    patches[1].Transpose(); //10
    patches[2].Transpose(); //01
    patches[3].Transpose(); //11
}

void BezierPatch::SplitU(BezierPatch patches[2], float u)
{
    for (int i = 0; i < N; ++i) {
        bezierSplit(&patches[0].m_cp[i*N + 0],
            &patches[1].m_cp[i*N + 0],
            &m_cp[i*N + 0], u);
    }
}

void BezierPatch::SplitV(BezierPatch patches[2], float v)
{
    for (int i = 0; i < N; ++i) {
        bezierSplitV(&patches[0].m_cp[i],
            &patches[1].m_cp[i],
            &m_cp[i], v);
    }
}

void BezierPatch::CropU(BezierPatch patch, float u0, float u1)
{
    for (int i = 0; i < N; ++i) {
        bezierCrop(&patch.m_cp[i*N + 0], &m_cp[i*N + 0], u0, u1);
    }
}

void BezierPatch::CropV(BezierPatch patch, float v0, float v1)
{
    for (int i = 0; i < N; ++i) {
        bezierCropV(&patch.m_cp[i], &m_cp[i], v0, v1);
    }
}

bool BezierPatch::Crop(BezierPatch patch, float u0, float u1, float v0, float v1)
{
    float3 tmp[Ncp];
    for (int i = 0; i < N; ++i) bezierCrop(&tmp[i*N + 0], &m_cp[i*N + 0], u0, u1);
    for (int i = 0; i < N; ++i) bezierCropV(&patch.m_cp[i], &tmp[i], v0, v1);
    return true;
}

float3 BezierPatch::GetLv()
{
    return Get(0, N - 1) - Get(0, 0) + Get(N - 1, N - 1) - Get(N - 1, 0);
}

float3 BezierPatch::GetLu()
{
    return Get(N - 1, 0) - Get(0, 0) + Get(N - 1, N - 1) - Get(0, N - 1);
}
*/

#endif // IstBezierPatch_h
