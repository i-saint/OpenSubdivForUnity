#include "pch.h"
#include "BezierPatch.h"



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
    for (int i = 0; i < 16; ++i)
    {
        m_cp[i] = p[i];
    }
}

BezierPatch::BezierPatch(const BezierPatch &v)
{
    for (int i = 0; i < 16; ++i)
    {
        m_cp[i] = v[i];
    }
}

BezierPatch::BezierPatch(const BezierPatch &v, const float3x3 &m)
{
    for (int i = 0; i < 16; ++i)
    {
        m_cp[i] = m * v[i];
    }
}

BezierPatch::BezierPatch(const BezierPatch &v, const float4x4 &m)
{
    for (int i = 0; i < 16; ++i)
    {
        m_cp[i] = float3(m * float4(v[i], 1.0f));
    }
}

BezierPatch::BezierPatch(const BezierPatch &v, float u0, float u1, float v0, float v1)
{
    for (int i = 0; i < 16; ++i) { m_cp[i] = v[i]; }

    BezierPatch tmp;
    CropU(tmp, u0, u1);
    tmp.CropV(*this, v0, v1);
}

const float3& BezierPatch::Get(int i, int j) const
{
    return m_cp[4*j + i];
}

void BezierPatch::Set(int i, int j, const float3& v)
{
    m_cp[4*j + i] = v;
}


float3 BezierPatch::GetMin() const
{
    float3 min_ = m_cp[0];
    for (int i = 1; i < 16; ++i) {
        min_ = min(min_, m_cp[i]);
    }
    return min_;
}

float3 BezierPatch::GetMax() const
{
    float3 max_ = m_cp[0];
    for (int i = 1; i < 16; ++i) {
        max_ = max(max_, m_cp[i]);
    }
    return max_;
}

void BezierPatch::GetMinMax(float3 &o_min, float3 &o_max, float epsilon) const
{
    o_min = o_max = m_cp[0];
    for (int i = 1; i < 16; ++i)
    {
        o_min = min(o_min, m_cp[i]);
        o_max = max(o_max, m_cp[i]);
    }
    o_min -= epsilon;
    o_max += epsilon;
}


float3 BezierPatch::GetRoughNormal() const
{
    float3 LU = m_cp[ 3] - m_cp[0];
    float3 LV = m_cp[12] - m_cp[0];
    return normalize(cross(LV, LU));
}

float3 BezierPatch::Evaluate(const float2& uv) const
{
    float3 b[4];
    for (int i = 0; i < 4; ++i) {
        b[i] = BPEvaluate(uv.x, m_cp + i*4);
    }
    return BPEvaluate(uv.y, b);
}

float3 BezierPatch::EvaluateDu(const float2& uv) const
{
    float3 b[4];
    for (int i = 0; i < 4; ++i) {
        b[i] = BPEvaluateD(uv.x, m_cp + i*4);
    }
    return BPEvaluate(uv.y, b);
}

float3 BezierPatch::EvaluateDv(const float2& uv) const
{
    float3 b[4];
    for (int i = 0; i < 4; ++i) {
        b[i] = BPEvaluate(uv.x, m_cp + i*4);
    }
    return BPEvaluateD(uv.y, b);
}

float3 BezierPatch::EvaluateNormal(const float2& uv) const
{
    float3 du = EvaluateDu(uv);
    float3 dv = EvaluateDv(uv);
    return normalize(cross(dv, du));
}

void BezierPatch::Split(BezierPatch &dst0, BezierPatch &dst1, BezierPatch &dst2, BezierPatch &dst3, const float2& uv) const
{
    BezierPatch tmp0, tmp1;

    // split U
    SplitU(tmp0, tmp1, uv.x);

    // uv -> vu
    tmp0.Transpose(); // 00 01
    tmp1.Transpose(); // 10 11

    // split V
    tmp0.SplitU(dst0, dst2, uv.y);
    tmp1.SplitU(dst1, dst3, uv.y);

    // vu -> uv
    dst0.Transpose(); //00
    dst1.Transpose(); //10
    dst2.Transpose(); //01
    dst3.Transpose(); //11
}

void BezierPatch::SplitU(BezierPatch &dst0, BezierPatch &dst1, float u) const
{
    float t = u;
    const float3 *src = m_cp;
    for (int i = 0; i < 16; i += 4) {
        float S = 1.0f - t;
        float3 p0 = src[i + 0];
        float3 p1 = src[i + 1];
        float3 p2 = src[i + 2];
        float3 p3 = src[i + 3];
        dst0[i + 0] = p0;
        dst0[i + 1] = p0*S + p1*t;
        dst0[i + 2] = p0*S*S + p1 * 2.0f * S*t + p2*t*t;
        dst0[i + 3] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;
        dst1[i + 0] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;
        dst1[i + 1] = p3*t*t + p2 * 2.0f * t*S + p1*S*S;
        dst1[i + 2] = p3*t + p2*S;
        dst1[i + 3] = p3;
    }
}

void BezierPatch::SplitV(BezierPatch &dst0, BezierPatch &dst1, float v) const
{
    float t = v;
    const float3 *src = m_cp;
    for (int i = 0; i < 4; ++i) {
        float S = 1.0f - t;
        float3 p0 = src[i + 0];
        float3 p1 = src[i + 4];
        float3 p2 = src[i + 8];
        float3 p3 = src[i +12];
        dst0[i + 0] = p0;
        dst0[i + 4] = p0*S + p1*t;
        dst0[i + 8] = p0*S*S + p1 * 2.0f * S*t + p2*t*t;
        dst0[i +12] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;
        dst1[i + 0] = p0*S*S*S + p1 * 3.0f * S*S*t + p2 * 3.0f * S*t*t + p3*t*t*t;
        dst1[i + 4] = p3*t*t + p2 * 2.0f * t*S + p1*S*S;
        dst1[i + 8] = p3*t + p2*S;
        dst1[i +12] = p3;
    }
}

void BezierPatch::Crop(BezierPatch &dst, const float2& uv0, const float2& uv1) const
{
    BezierPatch tmp;
    CropU(tmp, uv0.x, uv1.x);
    tmp.CropV(dst, uv0.y, uv1.y);
}

void BezierPatch::CropU(BezierPatch &dst, float u0, float u1) const
{
    float s = u0;
    float t = u1;
    const float3 *src = m_cp;
    for (int i = 0; i < 16; i += 4) {
        const float3 p0 = src[i + 0];
        const float3 p1 = src[i + 1];
        const float3 p2 = src[i + 2];
        const float3 p3 = src[i + 3];
        float T = 1.0f - s;
        float S = 1.0f - t;
        s = 1.0f - T;
        t = 1.0f - S;
        dst[i + 0] = (p0*(T*T)*T + p3*(s*s)*s) + (p1*(s*T)*(3.0f * T) + p2*(s*s)*(3.0f * T));
        dst[i + 1] = (p0*(T*T)*S + p3*(s*s)*t) + (p1*T*(2.0f * (S*s) + T*t) + p2*s*(2.0f * (t*T) + (s*S)));
        dst[i + 2] = (p3*(t*t)*s + p0*(S*S)*T) + (p2*t*(2.0f * (s*S) + t*T) + p1*S*(2.0f * (T*t) + (S*s)));
        dst[i + 3] = (p3*(t*t)*t + p0*(S*S)*S) + (p2*(S*t)*(3.0f * t) + p1*(S*S)*(3.0f * t));
    }
}

void BezierPatch::CropV(BezierPatch &dst, float v0, float v1) const
{
    float s = v0;
    float t = v1;
    const float3 *src = m_cp;
    for (int i = 0; i < 4; ++i) {
        const float3 p0 = src[i + 0];
        const float3 p1 = src[i + 4];
        const float3 p2 = src[i + 8];
        const float3 p3 = src[i +12];
        float T = 1.0f - s;
        float S = 1.0f - t;
        s = 1.0f - T;
        t = 1.0f - S;
        dst[i + 0] = (p0*(T*T)*T + p3*(s*s)*s) + (p1*(s*T)*(3.0f * T) + p2*(s*s)*(3.0f * T));
        dst[i + 4] = (p0*(T*T)*S + p3*(s*s)*t) + (p1*T*(2.0f * (S*s) + T*t) + p2*s*(2.0f * (t*T) + (s*S)));
        dst[i + 8] = (p3*(t*t)*s + p0*(S*S)*T) + (p2*t*(2.0f * (s*S) + t*T) + p1*S*(2.0f * (T*t) + (S*s)));
        dst[i +12] = (p3*(t*t)*t + p0*(S*S)*S) + (p2*(S*t)*(3.0f * t) + p1*(S*S)*(3.0f * t));
    }
}

float3 BezierPatch::GetLv() const
{
    return Get(0, 4 - 1) - Get(0, 0) + Get(4 - 1, 4 - 1) - Get(4 - 1, 0);
}

float3 BezierPatch::GetLu() const
{
    return Get(4 - 1, 0) - Get(0, 0) + Get(4 - 1, 4 - 1) - Get(0, 4 - 1);
}



void BezierPatch::Transpose()
{
    std::swap(m_cp[1 * 4 + 0], m_cp[0 * 4 + 1]);
    std::swap(m_cp[2 * 4 + 0], m_cp[0 * 4 + 2]);
    std::swap(m_cp[3 * 4 + 0], m_cp[0 * 4 + 3]);
    std::swap(m_cp[2 * 4 + 1], m_cp[1 * 4 + 2]);
    std::swap(m_cp[3 * 4 + 1], m_cp[1 * 4 + 3]);
    std::swap(m_cp[3 * 4 + 2], m_cp[2 * 4 + 3]);
}

void BezierPatch::Transform(const float3x3 &m)
{
    ForEach([&m](float3 &cp) { cp = m * cp; });
}

void BezierPatch::Transform(const float4x4 &m)
{
    ForEach([&m](float3 &cp) { cp = float3(m * float4(cp, 1.0f)); });
}
