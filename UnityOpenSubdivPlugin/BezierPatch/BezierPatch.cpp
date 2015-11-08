#include "pch.h"
#include "BezierPatch.h"



static inline void BPSplitU(float3 a[4], float3 b[4], const float3 p[4], float t)
{
    float S = 1.0f - t;
    float3 p0 = p[0];
    float3 p1 = p[1];
    float3 p2 = p[2];
    float3 p3 = p[3];
    a[0] = p0;
    a[1] = p0*S + p1*t;
    a[2] = p0*S*S + p1 * 2.0f * S*t + p2*t*t;
    a[3] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;

    b[0] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;
    b[1] = p3*t*t + p2 * 2.0f * t*S + p1*S*S;
    b[2] = p3*t + p2*S;
    b[3] = p3;
}

static inline void BPSplitV(float3 a[16], float3 b[16], const float3 p[16], float t)
{
    float S = 1.0f - t;
    float3 p0 = p[ 0];
    float3 p1 = p[ 4];
    float3 p2 = p[ 8];
    float3 p3 = p[12];
    a[ 0] = p0;
    a[ 4] = p0*S + p1*t;
    a[ 8] = p0*S*S + p1 * 2.0f * S*t + p2*t*t;
    a[12] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;

    b[ 0] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;
    b[ 4] = p3*t*t + p2 * 2.0f * t*S + p1*S*S;
    b[ 8] = p3*t + p2*S;
    b[12] = p3;
}

static inline void BPCropU(float3 out[4], const float3 in[4], float s, float t)
{
    const float3 p0 = in[0];
    const float3 p1 = in[1];
    const float3 p2 = in[2];
    const float3 p3 = in[3];
    float T = 1.0f - s;
    float S = 1.0f - t;
    s = 1.0f - T;
    t = 1.0f - S;
    out[0] = (p0*(T*T)*T + p3*(s*s)*s) + (p1*(s*T)*(3.0f * T) + p2*(s*s)*(3.0f * T));
    out[1] = (p0*(T*T)*S + p3*(s*s)*t) + (p1*T*(2.0f * (S*s) + T*t) + p2*s*(2.0f * (t*T) + (s*S)));
    out[2] = (p3*(t*t)*s + p0*(S*S)*T) + (p2*t*(2.0f * (s*S) + t*T) + p1*S*(2.0f * (T*t) + (S*s)));
    out[3] = (p3*(t*t)*t + p0*(S*S)*S) + (p2*(S*t)*(3.0f * t) + p1*(S*S)*(3.0f * t));
}

static inline void BPCropV(float3 out[16], const float3 in[16], float s, float t)
{
    const float3 p0 = in[ 0];
    const float3 p1 = in[ 4];
    const float3 p2 = in[ 8];
    const float3 p3 = in[12];
    float T = 1.0f - s;
    float S = 1.0f - t;
    s = 1.0f - T;
    t = 1.0f - S;
    out[ 0] = (p0*(T*T)*T + p3*(s*s)*s) + (p1*(s*T)*(3 * T) + p2*(s*s)*(3.0f * T));
    out[ 4] = (p0*(T*T)*S + p3*(s*s)*t) + (p1*T*(2.0f * (S*s) + T*t) + p2*s*(2.0f * (t*T) + (s*S)));
    out[ 8] = (p3*(t*t)*s + p0*(S*S)*T) + (p2*t*(2.0f * (s*S) + t*T) + p1*S*(2.0f * (T*t) + (S*s)));
    out[12] = (p3*(t*t)*t + p0*(S*S)*S) + (p2*(S*t)*(3.0f * t) + p1*(S*S)*(3.0f * t));
}

static inline float3 BPEvaluate(float t, const float3 cp[4])
{
    float it = 1.0f - t;
    return cp[0] * (it*it*it)
         + cp[1] * (3.0f*(it*it*t))
         + cp[2] * (3.0f*(it*t*t))
         + cp[3] * (t*t*t);
}

static inline float3 BPEvaluateD(float t, const float3 cp[4])
{
    float t2 = t * t;
    return cp[0] * (3.0f * t2 *-1.0f + 2.0f * t * 3.0f - 3.0f)
         + cp[1] * (3.0f * t2 * 3.0f + 2.0f * t *-6.0f + 3.0f)
         + cp[2] * (3.0f * t2 *-3.0f + 2.0f * t * 3.0f )
         + cp[3] * (3.0f * t2 * 1.0f );
}



BezierPatch::BezierPatch()
{
}

BezierPatch::BezierPatch(const float3 *p)
{
    for (int i = 0; i < Ncp; ++i)
    {
        m_cp[i] = p[i];
    }
}

BezierPatch::BezierPatch(const BezierPatch &v)
{
    for (int i = 0; i < Ncp; ++i)
    {
        m_cp[i] = v[i];
    }
}

BezierPatch::BezierPatch(const BezierPatch &v, const float3x3 &m)
{
    for (int i = 0; i < Ncp; ++i)
    {
        m_cp[i] = m * v[i];
    }
}

BezierPatch::BezierPatch(const BezierPatch &v, const float4x4 &m)
{
    for (int i = 0; i < Ncp; ++i)
    {
        m_cp[i] = float3(m * float4(v[i], 1.0f));
    }
}

BezierPatch::BezierPatch(const BezierPatch &v, float u0, float u1, float v0, float v1)
{
    for (int i = 0; i < Ncp; ++i) { m_cp[i] = v[i]; }

    float3 tmp[Ncp];
    for (int i = 0; i < N; ++i) { BPCropU(&tmp[i*N + 0], &m_cp[i*N + 0], u0, u1); }
    for (int i = 0; i < N; ++i) { BPCropV(&m_cp[i], &tmp[i], v0, v1); }
}

const float3& BezierPatch::Get(int i, int j) const
{
    return m_cp[N*j + i];
}

void BezierPatch::Set(int i, int j, const float3& v)
{
    m_cp[N*j + i] = v;
}


float3 BezierPatch::GetMin() const
{
    float3 min_ = m_cp[0];
    for (int i = 1; i < Ncp; ++i) {
        min_ = min(min_, m_cp[i]);
    }
    return min_;
}

float3 BezierPatch::GetMax() const
{
    float3 max_ = m_cp[0];
    for (int i = 1; i < Ncp; ++i) {
        max_ = max(max_, m_cp[i]);
    }
    return max_;
}

void BezierPatch::GetMinMax(float3 &o_min, float3 &o_max, float epsilon) const
{
    o_min = o_max = m_cp[0];
    for (int i = 1; i < Ncp; ++i)
    {
        o_min = min(o_min, m_cp[i]);
        o_max = max(o_max, m_cp[i]);
    }
    o_min -= epsilon;
    o_max += epsilon;
}


float3 BezierPatch::Evaluate(const float2& uv) const
{
    float3 b[N];
    for (int i = 0; i < N; ++i) {
        b[i] = BPEvaluate(uv.x, m_cp + i*N);
    }
    return BPEvaluate(uv.y, b);
}

float3 BezierPatch::EvaluateDu(const float2& uv) const
{
    float3 b[N];
    for (int i = 0; i < N; ++i) {
        b[i] = BPEvaluateD(uv.x, m_cp + i*N);
    }
    return BPEvaluate(uv.y, b);
}

float3 BezierPatch::EvaluateDv(const float2& uv) const
{
    float3 b[N];
    for (int i = 0; i < N; ++i) {
        b[i] = BPEvaluate(uv.x, m_cp + i*N);
    }
    return BPEvaluateD(uv.y, b);
}

float3 BezierPatch::EvaluateNormal(const float2& uv) const
{
    float3 du = EvaluateDu(uv);
    float3 dv = EvaluateDv(uv);
    return normalize(cross(dv, du));
}

void BezierPatch::Split(BezierPatch patches[4], float u, float v) const
{
    BezierPatch tmp0, tmp1;

    // split U
    for (int i = 0; i < N; ++i) {
        BPSplitU(&tmp0[i*N + 0], &tmp1[i*N + 0], &m_cp[i*N + 0], u);
    }

    // uv -> vu
    tmp0.Transpose(); // 00 01
    tmp1.Transpose(); // 10 11

    // split V
    for (int i = 0; i < N; ++i) {
        BPSplitU(&patches[0].m_cp[i*N + 0],
            &patches[2].m_cp[i*N + 0],
            &tmp0[i*N + 0], v);
        BPSplitU(&patches[1].m_cp[i*N + 0],
            &patches[3].m_cp[i*N + 0],
            &tmp1[i*N + 0], v);
    }

    // vu -> uv
    patches[0].Transpose(); //00
    patches[1].Transpose(); //10
    patches[2].Transpose(); //01
    patches[3].Transpose(); //11
}

void BezierPatch::SplitU(BezierPatch patches[2], float u) const
{
    for (int i = 0; i < N; ++i) {
        BPSplitU(&patches[0].m_cp[i*N + 0],
            &patches[1].m_cp[i*N + 0],
            &m_cp[i*N + 0], u);
    }
}

void BezierPatch::SplitV(BezierPatch patches[2], float v) const
{
    for (int i = 0; i < N; ++i) {
        BPSplitV(&patches[0].m_cp[i],
            &patches[1].m_cp[i],
            &m_cp[i], v);
    }
}

void BezierPatch::CropU(BezierPatch &patch, float u0, float u1) const
{
    for (int i = 0; i < N; ++i) {
        BPCropU(&patch.m_cp[i*N + 0], &m_cp[i*N + 0], u0, u1);
    }
}

void BezierPatch::CropV(BezierPatch &patch, float v0, float v1) const
{
    for (int i = 0; i < N; ++i) {
        BPCropV(&patch.m_cp[i], &m_cp[i], v0, v1);
    }
}

bool BezierPatch::Crop(BezierPatch &patch, float u0, float u1, float v0, float v1) const
{
    float3 tmp[Ncp];
    for (int i = 0; i < N; ++i) BPCropU(&tmp[i*N + 0], &m_cp[i*N + 0], u0, u1);
    for (int i = 0; i < N; ++i) BPCropV(&patch.m_cp[i], &tmp[i], v0, v1);
    return true;
}

float3 BezierPatch::GetLv() const
{
    return Get(0, N - 1) - Get(0, 0) + Get(N - 1, N - 1) - Get(N - 1, 0);
}

float3 BezierPatch::GetLu() const
{
    return Get(N - 1, 0) - Get(0, 0) + Get(N - 1, N - 1) - Get(0, N - 1);
}



void BezierPatch::Transpose()
{
    std::swap(m_cp[1 * N + 0], m_cp[0 * N + 1]);
    std::swap(m_cp[2 * N + 0], m_cp[0 * N + 2]);
    std::swap(m_cp[3 * N + 0], m_cp[0 * N + 3]);
    std::swap(m_cp[2 * N + 1], m_cp[1 * N + 2]);
    std::swap(m_cp[3 * N + 1], m_cp[1 * N + 3]);
    std::swap(m_cp[3 * N + 2], m_cp[2 * N + 3]);
}

void BezierPatch::Transform(const float3x3 &m)
{
    ForEach([&m](float3 &cp) { cp = m * cp; });
}

void BezierPatch::Transform(const float4x4 &m)
{
    ForEach([&m](float3 &cp) { cp = float3(m * float4(cp, 1.0f)); });
}
