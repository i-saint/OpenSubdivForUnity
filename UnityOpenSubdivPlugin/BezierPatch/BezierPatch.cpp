#include "pch.h"
#include "BezierPatch.h"








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
    for (int i = 0; i < Ncp; ++i) m_cp[i] = v[i];
    float3 tmp[Ncp];
    for (int i = 0; i < N; ++i) bezierCrop(&tmp[i*N + 0], &m_cp[i*N + 0], u0, u1);
    for (int i = 0; i < N; ++i) bezierCropV(&m_cp[i], &tmp[i], v0, v1);
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
        b[i] = evaluate(uv.x, m_cp + i*N);
    }
    return evaluate(uv.y, b);
}

float3 BezierPatch::EvaluateDu(const float2& uv) const
{
    float3 b[N];
    for (int i = 0; i < N; ++i) {
        b[i] = evaluateD(uv.x, m_cp + i*N);
    }
    return evaluate(uv.y, b);
}

float3 BezierPatch::EvaluateDv(const float2& uv) const
{
    float3 b[N];
    for (int i = 0; i < N; ++i) {
        b[i] = evaluate(uv.x, m_cp + i*N);
    }
    return evaluateD(uv.y, b);
}

float3 BezierPatch::EvaluateNormal(const float2& uv) const
{
    float3 du = EvaluateDu(uv);
    float3 dv = EvaluateDv(uv);
    return normalize(cross(du, dv));
}

void BezierPatch::Split(BezierPatch patches[4], float u, float v) const
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

void BezierPatch::SplitU(BezierPatch patches[2], float u) const
{
    for (int i = 0; i < N; ++i) {
        bezierSplit(&patches[0].m_cp[i*N + 0],
            &patches[1].m_cp[i*N + 0],
            &m_cp[i*N + 0], u);
    }
}

void BezierPatch::SplitV(BezierPatch patches[2], float v) const
{
    for (int i = 0; i < N; ++i) {
        bezierSplitV(&patches[0].m_cp[i],
            &patches[1].m_cp[i],
            &m_cp[i], v);
    }
}

void BezierPatch::CropU(BezierPatch &patch, float u0, float u1) const
{
    for (int i = 0; i < N; ++i) {
        bezierCrop(&patch.m_cp[i*N + 0], &m_cp[i*N + 0], u0, u1);
    }
}

void BezierPatch::CropV(BezierPatch &patch, float v0, float v1) const
{
    for (int i = 0; i < N; ++i) {
        bezierCropV(&patch.m_cp[i], &m_cp[i], v0, v1);
    }
}

bool BezierPatch::Crop(BezierPatch &patch, float u0, float u1, float v0, float v1) const
{
    float3 tmp[Ncp];
    for (int i = 0; i < N; ++i) bezierCrop(&tmp[i*N + 0], &m_cp[i*N + 0], u0, u1);
    for (int i = 0; i < N; ++i) bezierCropV(&patch.m_cp[i], &tmp[i], v0, v1);
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




template<int N, int k>
struct Binomial
{
    const static int value = (Binomial<N - 1, k - 1>::value + Binomial<N - 1, k>::value);
};

template<>
struct Binomial<0, 0>
{
    const static int value = 1;
};
template<int N>
struct Binomial<N, 0>
{
    const static int value = 1;
};
template<int N>
struct Binomial<N, N>
{
    const static int value = 1;
};

template<int N, int k>
inline int binomial()
{
    return Binomial<N, k>::value;
}


template<int N, int I>
inline float bernstein(float t)
{
    return float(binomial<N, I>() * std::pow(1.0f - t, N - I) * std::pow(t, I));
}

template<int N, int I>
struct BezierCurveEvaluate
{
    static float3 evaluate(const float t, const float3 *const cp)
    {
        return cp[I] * bernstein<N, I>(t) +
            BezierCurveEvaluate<N, I - 1>::evaluate(t, cp);
    }
};

template<int N>
struct BezierCurveEvaluate<N, 0>
{
    static float3 evaluate(const float t, const float3 *const cp)
    {
        return cp[0] * bernstein<N, 0>(t);
    }
};



template <int N, int I, int STRIDE>
struct BezierSplit
{
    BezierSplit(float3 a[], float3 b[], const float3 p[], float t)
    {
        int tn = N - 1 - I;
        a[I*STRIDE] = p[0];
        b[tn*STRIDE] = p[tn*STRIDE];
        float3 tmp[(N - 1)*STRIDE];
        float s = 1 - t;
        t = 1 - s;
        for (int j = 0; j < tn; j++) {
            tmp[j*STRIDE] = s*p[j*STRIDE] + t*p[(j + 1)*STRIDE];
        }
        BezierSplit<N, I + 1, STRIDE>(a, b, &tmp[0], t);
    }
};

template <int N, int STRIDE>
struct BezierSplit<N, N, STRIDE>
{
    BezierSplit(float3 a[], float3 b[], const float3 p[], float t)
    {
    }
};



template <int N, int STRIDE>
struct BezierCrop
{
    BezierCrop(float3 r[], const float3 p[], float s, float t)
    {
        float3 tmp0[N*STRIDE], tmp1[N*STRIDE];
        BezierSplit<N, 0, STRIDE> split1(&tmp0[0], &tmp1[0], p, t);
        BezierSplit<N, 0, STRIDE> split2(&tmp1[0], r, &tmp0[0], s / t);
    }
};

template <int STRIDE>
struct BezierCrop<4, STRIDE>
{
    BezierCrop(float3 r[], const float3 p[], float s, float t)
    {
        const float3 &p0 = p[0 * STRIDE];
        const float3 &p1 = p[1 * STRIDE];
        const float3 &p2 = p[2 * STRIDE];
        const float3 &p3 = p[3 * STRIDE];
        // watertight cropping for cubic b-spline
        //
        //   p0    p1   p2    p3
        //   0----s--------t---1
        //   1----T--------S---0
        //
        float T = 1 - s;
        float S = 1 - t;
        s = 1 - T;
        t = 1 - S;
        r[0 * STRIDE] = (p0*(T*T)*T + p3*(s*s)*s) + (p1*(s*T)*(3 * T) + p2*(s*s)*(3 * T));
        r[1 * STRIDE] = (p0*(T*T)*S + p3*(s*s)*t) + (p1*T*(2 * (S*s) + T*t) + p2*s*(2 * (t*T) + (s*S)));
        r[2 * STRIDE] = (p3*(t*t)*s + p0*(S*S)*T) + (p2*t*(2 * (s*S) + t*T) + p1*S*(2 * (T*t) + (S*s)));
        r[3 * STRIDE] = (p3*(t*t)*t + p0*(S*S)*S) + (p2*(S*t)*(3 * t) + p1*(S*S)*(3 * t));
    }
};

void BezierPatch::bezierSplit(float3 a[], float3 b[], const float3 p[], float t) const
{
    BezierSplit<N, 0, 1> split(a, b, p, t);
}

void BezierPatch::bezierSplitV(float3 a[], float3 b[], const float3 p[], float t) const
{
    BezierSplit<N, 0, N> split(a, b, p, t);
}

void BezierPatch::bezierCrop(float3 r[], const float3 p[], float s, float t) const
{
    BezierCrop<N, 1> crop(r, p, s, t);
}

void BezierPatch::bezierCropV(float3 r[], const float3 p[], float s, float t) const
{
    BezierCrop<N, N> crop(r, p, s, t);
}


float3 BezierPatch::evaluate(float t, const float3 *cp)
{
    //return BezierCurveEvaluate<N - 1, N - 1>::evaluate(t, cp);

    float it = 1.0f - t;
    return cp[0] * (       (it*it*it)          )
         + cp[1] * (3.0f * (it*it   ) * (t)    )
         + cp[2] * (3.0f * (it      ) * (t*t)  )
         + cp[3] * (                    (t*t*t));
}

float3 BezierPatch::evaluateD(float t, const float3 *cp)
{
    float t2 = t * t;
    return cp[0] * (3.0f * t2 *-1.0f + 2.0f * t * 3.0f - 3.0f)
         + cp[1] * (3.0f * t2 * 3.0f + 2.0f * t *-6.0f + 3.0f)
         + cp[2] * (3.0f * t2 *-3.0f + 2.0f * t * 3.0f       )
         + cp[3] * (3.0f * t2 * 1.0f                         );
}




void BezierPatch::Transpose()
{
    for (int y = 0; y < N; ++y) {
        for (int x = y + 1; x < N; ++x) {
            std::swap(m_cp[x*N + y], m_cp[y*N + x]);
        }
    }
}

void BezierPatch::Transform(const float3x3 &m)
{
    ForEach([&m](float3 &cp) { cp = m * cp; });
}

void BezierPatch::Transform(const float4x4 &m)
{
    ForEach([&m](float3 &cp) { cp = float3(m * float4(cp, 1.0f)); });
}
