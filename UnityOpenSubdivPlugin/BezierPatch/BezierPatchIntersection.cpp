#include "pch.h"
#include "BezierPatch.h"
#include "BezierPatchIntersection.h"

#define trace(...)
#define N 4


struct Epsilon
{
    static float GetEps() { return 1e-4f; }
    static float GetUvEps() { return 1.0f / 32.0f; }
};


float3x3 Rotate2D(const float3 &dx)
{
    float2 x = normalize(float2(dx.x, dx.y));
    float2 y = float2(-x[1], x[0]);
    return float3x3(
        x[0], x[1], 0.0f,
        y[0], y[1], 0.0f,
        0.0f, 0.0f, 1.0f );
}


float4x4 RayTransform(const float3& pos, const float3& dir)
{
    float3 z = dir;
    int plane = 0;
    if (fabs(z[1]) < fabs(z[plane])) plane = 1;
    if (fabs(z[2]) < fabs(z[plane])) plane = 2;

    float3 x = (plane == 0) ? float3(1.0f, 0.0f, 0.0f) :
        ((plane == 1) ? float3(0.0f, 1.0f, 0.0f) :
            float3(0.0f, 0.0f, 1.0f));
    float3 y = normalize(cross(z, x));
    x = cross(y, z);

    float4x4 rot = float4x4(
        x[0], x[1], x[2], 0.0f,
        y[0], y[1], y[2], 0.0f,
        z[0], z[1], z[2], 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f);
    float4x4 trs = float4x4(
        1.0f, 0.0f, 0.0f, -pos[0],
        0.0f, 1.0f, 0.0f, -pos[1],
        0.0f, 0.0f, 1.0f, -pos[2],
        0.0f, 0.0f, 0.0f, 1.0f);
    return rot*trs;
}



void RotateU(BezierPatch& patch)
{
    float3 dx = patch.GetLv();
    patch.Transform(Rotate2D(dx));
}

void RotateV(BezierPatch& patch)
{
    float3 dx = patch.GetLu();
    patch.Transform(Rotate2D(dx));
}

bool IsEps(float3 const & min, float3 const & max, float eps)
{
    //return ((max[1]-min[1])<=eps);
    float xw = max[0] - min[0];
    float yw = max[1] - min[1];
    return (xw <= eps && yw <= eps);
}

bool IsLevel(int level, int max_level)
{
    return (level >= max_level);
}

bool IsClip(int level)
{
    static const int div_level[] =
    {
        1,
        1,1,
        2,2,//4
        3,3,3,3,//8
        4,4,4,4,4,4,4,4,//16
    };
    int nlevel = 0;
    if (N <= 16) {
        nlevel = div_level[N];
    }
    else {
        nlevel = (int)ceil(log2(N));
    }
    return nlevel * 2 <= level;
}

float Lerp(float a, float b, float u)
{
    float t = float(1) - (float(1) - u);
    float s = float(1) - t;
    return s*a + t*b;
}

static void CoarseSort(int order[2], BezierPatch tmp[2])
{
#if USE_COARSESORT
    float zs[2];
    for (int i = 0; i < 2; ++i) {
        int u = BezierPatch::Order >> 1;
        int v = BezierPatch::Order >> 1;

        zs[i] = tmp[i].Get(u, v)[2];
    }
    if (zs[0] <= zs[1]) {
        order[0] = 0;
        order[1] = 1;
    }
    else {
        order[0] = 1;
        order[1] = 0;
    }
#endif
}

// --------------------------------------------------------------------

struct Curve
{
    static const int LENGTH = 4;
    float2 operator[](int i) const { return v[i]; }
    float2 &operator[](int i) { return v[i]; }
    float2 v[LENGTH];
};

bool ScanMinMax(const Curve& p)
{
    int nUpper = 0, nLower = 0;
    for (int i = 0; i < Curve::LENGTH; ++i) {
        if (p[i][1] > 0) nUpper++;
        else nLower++;
        if (nUpper && nLower) return true;
    }
    return false;
}

template<typename T>
T Slope(T const a[], T const b[])
{
    T d0 = b[0] - a[0];
    T d1 = b[1] - a[1];
    return fabs(d0 / d1);
}

template<typename T>
T Dist(T const p0[], T const p1[])
{
    T a = fabs(p0[1]);
    T b = fabs(p1[1]);
    T h = a + b;
    T t = p0[0] + a*(p1[0] - p0[0]) / h;
    return t;
}

int ScanConvex(float ts[], const Curve & p)
{
    if (!ScanMinMax(p)) {
        return 0;
    }
    int n = 0;
    {
        int k = -1;
        int l = -1;
        float current = std::numeric_limits<float>::max();
        for (int i = 1; i < Curve::LENGTH; ++i) {
            if (p[i][1] * p[0][1] < 0) {
                float s = Slope(&p[i][0], &p[0][0]);
                if (s < current) {
                    current = s;
                    k = i;
                }
            }
        }
        if (k >= 0) {
            current = 0;
            for (int i = 0; i < k; ++i) {
                if (p[i][1] * p[k][1] < 0) {
                    float s = Slope(&p[i][0], &p[k][0]);
                    if (current < s) {
                        current = s;
                        l = i;
                    }
                }
            }
            if (l >= 0) {
                ts[n++] = Dist(&p[l][0], &p[k][0]);
            }
        }
    }
    {
        int k = -1;
        int l = -1;
        float current = std::numeric_limits<float>::max();
        for (int i = 0; i < Curve::LENGTH - 1; ++i) {
            if (p[i][1] * p[Curve::LENGTH - 1][1] < 0) {
                float s = Slope(&p[i][0], &p[Curve::LENGTH - 1][0]);
                if (s < current) {
                    current = s;
                    k = i;
                }
            }
        }
        if (k >= 0) {
            current = 0;
            for (int i = k + 1; i < Curve::LENGTH; ++i) {
                if (p[i][1] * p[k][1] < 0) {
                    float s = Slope(&p[i][0], &p[k][0]);
                    if (current < s) {
                        current = s;
                        l = i;
                    }
                }
            }
            if (l >= 0) {
                ts[n++] = Dist(&p[k][0], &p[l][0]);
            }
        }
    }
    return n;
}

bool XCheck(float rng[2], const Curve& curve)
{
    float t[4] = { 0 };
    int nn = ScanConvex(t, curve);
    if (nn) {
        float t0 = t[0];
        float t1 = t[0];

        for (int i = 1; i < nn; i++) {
            t0 = std::min(t0, t[i]);
            t1 = std::max(t1, t[i]);
        }

        rng[0] = t0;
        rng[1] = t1;
        return true;
    }
    return false;
}

bool GetRangeU(float out[2], BezierPatch const & patch)
{
    float t0 = 1;
    float t1 = 0;
    Curve curve;
    float delta = 1.0 / (N - 1);
    for (int i = 0; i < N; ++i) {
        curve[i][0] = i*delta;
    }
    float rng[2];
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            curve[i][1] = (patch.Get(i, j))[1];
        }
        if (XCheck(rng, curve)) {
            t0 = std::min(t0, rng[0]);
            t1 = std::max(t1, rng[1]);
        }
        else {
            return false;
        }
    }
    out[0] = t0;
    out[1] = t1;
    return true;
}

bool GetRangeV(float out[2], const BezierPatch& patch)
{
    float t0 = 1;
    float t1 = 0;
    Curve curve;
    float delta = 1.0 / (N - 1);
    for (int i = 0; i < N; ++i) {
        curve[i][0] = i*delta;
    }
    float rng[2];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            curve[j][1] = (patch.Get(i, j))[1];
        }
        if (XCheck(rng, curve)) {
            t0 = std::min(t0, rng[0]);
            t1 = std::max(t1, rng[1]);
        }
        else {
            return false;
        }
    }
    out[0] = t0;
    out[1] = t1;
    return true;
}


template<class T>
static bool SolveBilinearPatch(T* t, T* u, T* v, T tmin, T tmax, T tt, T uu, T vv) {
    if (tmin <= tt && tt <= tmax) {
        *t = tt;
        *u = uu;
        *v = vv;
        return true;
    }
    return false;
}

template <class T>
static int Solve2(T root[2], T const coeff[3]) {
    T A = coeff[0];
    T B = coeff[1];
    T C = coeff[2];
    if (fabs(A) <= std::numeric_limits<T>::min()) {
        if (fabs(B) <= std::numeric_limits<T>::min()) return 0;
        T x = -C / B;
        root[0] = x;
        return 1;
    }
    else {
        T D = B*B - 4 * A*C;
        if (D < 0) {
            return 0;
        }
        else if (D == 0) {
            T x = -0.5*B / A;
            root[0] = x;
            return 1;
        }
        else {
            T x1 = (fabs(B) + sqrt(D)) / (2.0 * A);
            if (B >= 0) {
                x1 = -x1;
            }
            T x2 = C / (A * x1);
            if (x1 > x2) std::swap(x1, x2);
            root[0] = x1;
            root[1] = x2;
            return 2;
        }
    }
}
template<class T>
static T ComputeU(T A1, T A2, T B1, T B2, T C1, T C2, T D1, T D2, T v) {
    //return div((v*(C1-C2)+(D1-D2)),(v*(A2-A1)+(B2-B1)));
    T a = v*A2 + B2;
    T b = v*(A2 - A1) + B2 - B1;
    if (fabs(b) >= fabs(a)) {
        return (v*(C1 - C2) + D1 - D2) / b;
    }
    else {
        return (-v*C2 - D2) / a;
    }
}
template<class T>
static T ComputeT(T a, T b, T c, T d, T iq, T u, T v)
{
    return ((u*v)*a + u*b + v*c + d)*iq;
}

template<typename T, typename V>
static bool TestBilinearPatch(T *t, T *u, T *v,
    V const p[4], T tmin, T tmax, float uvMargin)
{
    trace("Test bilinear (%f, %f, %f) (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n",
        p[0][0], p[0][1], p[0][2],
        p[1][0], p[1][1], p[1][2],
        p[2][0], p[2][1], p[2][2],
        p[3][0], p[3][1], p[3][2]);

    V const & p00 = p[0];
    V const & p10 = p[1];
    V const & p01 = p[2];
    V const & p11 = p[3];

    static const int nPlane = 2;
    V a = p11 - p10 - p01 + p00;
    V b = p10 - p00;
    V c = p01 - p00;
    V d = p00;

    //xz-zx
    float A1 = a[0];
    float B1 = b[0];
    float C1 = c[0];
    float D1 = d[0];

    //yz-zy
    float A2 = a[1];
    float B2 = b[1];
    float C2 = c[1];
    float D2 = d[1];

    float F1 = A2*C1 - A1*C2;
    float F2 = A2*D1 - A1*D2 + B2*C1 - B1*C2;
    float F3 = B2*D1 - B1*D2;

    float coeff[] = { F1,F2,F3 };
    float root[2] = {};
    //    printf("<%f, %f, %f>\n", F1, F2, F3);
    int nRet = Solve2(root, coeff);

    if (nRet) {
        bool bRet = false;
        trace("Solve bilinear:");
        for (int i = 0; i < nRet; ++i) {
            float vv = root[i];
            trace(" vv = %f, ", vv);
            if (0 - uvMargin <= vv && vv <= 1 + uvMargin) {//TODO
                vv = std::max(float(0), std::min(vv, float(1)));
                float uu = ComputeU(A1, A2, B1, B2, C1, C2, D1, D2, vv);
                trace(" uu = %f, ", uu);
                if (0 - uvMargin <= uu && uu <= 1 + uvMargin) {//TODO
                    uu = std::max(float(0), std::min(uu, float(1)));
                    float tt = ComputeT(a[nPlane], b[nPlane], c[nPlane], d[nPlane],
                        float(1), uu, vv);
                    trace("uv=(%f, %f)", uu, vv);
                    if (SolveBilinearPatch(t, u, v, tmin, tmax, tt, uu, vv)) {
                        tmax = *t;
                        bRet = true;
                    }
                }
            }
        }
        trace("\n");
        return bRet;
    }
    trace("\n");
    return false;
}

template <typename T, typename V>
inline T dot_(V const &a, V const &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template<typename T, typename V>
static bool TestTriangle(T *tt, T *uu, T *vv,
    V const & p0, V const & p1, V const & p2, V const & org, V const & dir, T tmin, T tmax)
{
    //typedef typename V::ElementType Scalar;
    T t, u, v;
    //-e1 = p0-p1
    V e1(p0 - p1);//vA
                  //-e2 = p0-p2
    V e2(p0 - p2);//vB 
                  //dir = GHI 
    V bDir(cross(e2, dir));

    T iM = dot_<T, V>(e1, bDir);

    if (0 < iM) {
        //p0-org
        V vOrg(p0 - org);

        u = dot_<T, V>(vOrg, bDir);
        if (u < 0 || iM < u)return false;

        V vE(cross(e1, vOrg));

        v = dot_<T, V>(dir, vE);
        if (v < 0 || iM < u + v)return false;

        t = -dot_<T, V>(e2, vE);
    }
    else if (iM < 0) {
        //p0-org
        V vOrg(p0 - org);//JKL

        u = dot_<T, V>(vOrg, bDir);
        if (u > 0 || iM > u)return false;

        V vE(cross(e1, vOrg));

        v = dot_<T, V>(dir, vE);
        if (v > 0 || iM > u + v)return false;

        t = -dot_<T, V>(e2, vE);
    }
    else {
        return false;
    }

    iM = T(1) / iM;
    t *= iM;
    if (t < tmin || tmax < t)return false;
    u *= iM;
    v *= iM;

    *tt = t;
    *uu = u;
    *vv = v;

    return true;
}


template<typename T, typename V>
static bool TestQuadPlane(T *t, T *u, T *v,
    V const p[4], V const & org, V const & dir, T tmin, T tmax)
{
    //typedef typename V::ElementType Scalar;

    V const & p00 = p[0];
    V const & p10 = p[1];
    V const & p01 = p[2];
    V const & p11 = p[3];

    T tt, uu, vv;
    bool bRet = false;
    if (TestTriangle(&tt, &uu, &vv, p00, p01, p10, org, dir, tmin, tmax))
    {
        T ww = T(1) - (uu + vv);
        *u = ww*T(0) + uu*T(0) + vv*T(1);//00 - 01 - 10
        *v = ww*T(0) + uu*T(1) + vv*T(0);//00 - 01 - 10
        *t = tt;
        tmax = tt;
        bRet = true;
    }
    if (TestTriangle(&tt, &uu, &vv, p10, p01, p11, org, dir, tmin, tmax))
    {
        T ww = T(1) - (uu + vv);
        *u = ww*T(1) + uu*T(0) + vv*T(1);//10 - 01 - 11
        *v = ww*T(0) + uu*T(1) + vv*T(1);//10 - 01 - 11
        *t = tt;
        tmax = tt;
        bRet = true;
    }

    return bRet;
}


BezierPatchIntersectionImpl::BezierPatchIntersectionImpl(const BezierPatch &patch) :
    m_patch(patch),
    m_eps(Epsilon::GetEps()),
    m_maxLevel(DEFAULT_MAX_LEVEL),
    m_uvMargin(true),
    m_cropUV(true), m_useBezierClip(true),
    m_useTriangle(false),
    m_directBilinear(false),
    m_wcpFlag(0)
{

    m_uRange[0] = m_vRange[0] = 0;
    m_uRange[1] = m_vRange[1] = 1;

    m_count = 0;

    patch.GetMinMax(m_min, m_max, 1e-3f);
}

BezierPatchIntersectionImpl::~BezierPatchIntersectionImpl()
{
}

void BezierPatchIntersectionImpl::SetEpsilon(float eps)
{
    m_eps = std::max(std::numeric_limits<float>::epsilon() * 1e+1f, eps);
}
void BezierPatchIntersectionImpl::SetMaxLevel(int maxLevel)
{
    m_maxLevel = maxLevel;
}
void BezierPatchIntersectionImpl::SetUVMergin(float uvMargin)
{
    m_uvMargin = uvMargin;
}
void BezierPatchIntersectionImpl::SetCropUV(bool cropUV)
{
    m_cropUV = cropUV;
}
void BezierPatchIntersectionImpl::SetUseBezierClip(bool useBezierClip)
{
    m_useBezierClip = useBezierClip;
}
void BezierPatchIntersectionImpl::SetUseTriangle(bool useTriangle)
{
    m_useTriangle = useTriangle;
}
void BezierPatchIntersectionImpl::SetDirectBilinear(bool directBilinear)
{
    m_directBilinear = directBilinear;
}
void BezierPatchIntersectionImpl::SetWatertightFlag(int wcpFlag)
{
    m_wcpFlag = wcpFlag;
}


bool BezierPatchIntersectionImpl::Test(BezierPatchHit &info, const Ray& r, float tmin, float tmax)
{
    return testInternal(info, r, tmin, tmax);
}


bool BezierPatchIntersectionImpl::testInternal(BezierPatchHit &info, const Ray& r, float tmin, float tmax)
{
    float4x4 mat = RayTransform(r.org, r.dir); // getZAlign

    bool bRet = false;
    UVT uvt;
    uvt.failFlag = 0;
    {
        BezierPatch patch(m_patch, mat);
        if (testBezierPatch(uvt, patch, tmin, tmax, m_eps)) {
            float t = uvt.t;
            float u = uvt.u;
            float v = uvt.v;

            u = m_uRange[0] * (1 - u) + m_uRange[1] * u;//global
            v = m_vRange[0] * (1 - v) + m_vRange[1] * v;//global
            info.t = t;
            info.u = u;
            info.v = v;
            info.clip_level = uvt.level;
            tmax = t;
            bRet = true;
            trace("hit t = %f, uv = (%f, %f)\n", t, u, v);
        }
    }

    if (uvt.failFlag == 0 or m_wcpFlag == 0) return bRet;

    {
        // TODO: still inefficient. wcpFlag knows which edge has to be split
        // and failFlag knows which edge is actually being tested.

        // test against split faces too
        /*
        +---+---+
        ^  | 2 | 3 |
        |  +---+---+
        |  | 0 | 1 |
        u  +---+---+
        v---->
        */
        static const float offsets[8] = {
            0.0f,0.0f,
            0.0f,0.5f,
            0.5f,0.0f,
            0.5f,0.5f
        };
        BezierPatch tmp[2];
        BezierPatch children[4];
        m_patch.SplitU(tmp, 0.5f);
        tmp[0].SplitV(&children[0], 0.5f);
        tmp[1].SplitV(&children[2], 0.5f);

        int ni = -1;
        for (int i = 0; i < 4; ++i) {
            children[i].Transform(mat);
            trace("Child pass %d\n", i);
            if (testBezierPatch(uvt, children[i], tmin, tmax, m_eps)) {
                float t = uvt.t;
                tmax = t;
                ni = i;
            }
        }

        if (0 <= ni) {
            int i = ni;
            float t = uvt.t;
            float u = uvt.u * 0.5f + offsets[2 * i + 0];
            float v = uvt.v * 0.5f + offsets[2 * i + 1];

            u = m_uRange[0] * (1.0f - u) + m_uRange[1] * u;//global
            v = m_vRange[0] * (1.0f - v) + m_vRange[1] * v;//global
            info.t = t;
            info.u = u;
            info.v = v;
            info.clip_level = uvt.level;
            bRet = true;
            trace("hit t = %f, uv = (%f, %f)\n", t, u, v);
        }
    }
    if (m_count > MAX_ITERATION) return false;
    return bRet;
}

bool BezierPatchIntersectionImpl::testBezierPatch(UVT& info, const BezierPatch& patch, float zmin, float zmax, float eps)
{
    float3 min, max;
    patch.GetMinMax(min, max, 1e-3f);

    if (0.0f < min.x || max.x < 0.0f) return false; //x
    if (0.0f < min.y || max.y < 0.0f) return false; //y
    if (max.z < zmin || zmax < min.z) return false; //z

    if (m_cropUV) return testBezierClipRangeU(info, patch, 0.0f, 1.0f, 0.0f, 1.0f, zmin, zmax, 0, m_maxLevel, eps);
    else          return testBezierClipU(info, patch, 0.0f, 1.0f, 0.0f, 1.0f, zmin, zmax, 0, m_maxLevel, eps);
}



bool BezierPatchIntersectionImpl::testBezierClipU(UVT& info, const BezierPatch & patch,
    float u0, float u1, float v0, float v1,
    float zmin, float zmax,
    int level, int max_level, float eps)
{
    BezierPatch tpatch(patch);
    RotateU(tpatch);

    if (++m_count > MAX_ITERATION || u0 > u1 || v0 > v1) {
        return false;
    }

    float3 min, max;
    patch.GetMinMax(min, max, 1e-3f);

    if (0.0f < min.x || max.x < 0.0f || // x
        0.0f < min.y || max.y < 0.0f || // y
        max.z < zmin || zmax < min.z )  // z
    {
#if USE_MINMAX_FAILFLAG
        // set failFlag only if it's very close
        if (fabs(min[0]) < eps ||
            fabs(max[0]) < eps ||
            fabs(min[1]) < eps ||
            fabs(max[1]) < eps) {
            int failFlag = ((u0 == float(0)) << 0)
                | ((u1 == float(1)) << 1)
                | ((v0 == float(0)) << 2)
                | ((v1 == float(1)) << 3);
            info->failFlag |= failFlag;
        }
#endif
        return false;
    }

    bool bClip = IsClip(level);
    if (bClip && (IsEps(min, max, eps) || IsLevel(level, max_level))) {
        return testBezierClipL(info, patch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1;
        float tt0 = 0;
        float tt1 = 1;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (GetRangeU(rng, tpatch)) {
                tt0 = rng[0];
                tt1 = rng[1];
                tw = tt1 - tt0;
            }
        }

        if (tw >= 0.4) {
            BezierPatch tmp[2];
            patch.SplitU(tmp, 0.5);
            float um = (u0 + u1)*0.5;
            float ut[] = { u0,um,um,u1 };

            int order[2] = { 0,1 };
            CoarseSort(order, tmp);
            bool bRet = false;
            if (ut[0] > ut[1]) std::swap(ut[0], ut[1]);
            if (ut[2] > ut[3]) std::swap(ut[2], ut[3]);
            for (int i = 0; i < 2; ++i) {
                if (testBezierClipV(info, tmp[order[i]], ut[2 * order[i]], ut[2 * order[i] + 1],
                    v0, v1, zmin, zmax, level + 1, max_level, eps)) {
                    zmax = info.t;
                    bRet = true;
                }
            }
            return bRet;
        }
        else {
            tt0 = std::max<float>(0.0, tt0 - Epsilon::GetUvEps());
            tt1 = std::min<float>(tt1 + Epsilon::GetUvEps(), 1.0);
            float ut[] = { Lerp(u0,u1,tt0),Lerp(u0,u1,tt1) };
            BezierPatch tmp;
            patch.CropU(tmp, tt0, tt1);
            return testBezierClipV(info, tmp, ut[0], ut[1],
                v0, v1, zmin, zmax, level + 1, max_level, eps);
        }
    }
    return false;
}

bool BezierPatchIntersectionImpl::testBezierClipV(UVT& info, const BezierPatch& patch,
    float u0, float u1, float v0, float v1, float zmin, float zmax,
    int level, int max_level, float eps)
{
    BezierPatch tpatch(patch);
    RotateV(tpatch);

    if (++m_count > MAX_ITERATION || u0 > u1 || v0 > v1) {
        return false;
    }

    trace("testBezierClipV (%f, %f) - (%f, %f) z:%f, %f  level=%d\n",
        u0, u1, v0, v1, zmin, zmax, level);

    float3 min, max;
    tpatch.GetMinMax(min, max, 1e-3f);

    if (0 < min[0] || max[0] < 0 || // x
        0 < min[1] || max[1] < 0 || // y
        max[2] < zmin || zmax < min[2]) { // z
#if USE_MINMAX_FAILFLAG
                                          // set failFlag only if it's very close
        if (fabs(min[0]) < eps ||
            fabs(max[0]) < eps ||
            fabs(min[1]) < eps ||
            fabs(max[1]) < eps) {
            int failFlag = ((u0 == float(0)) << 0)
                | ((u1 == float(1)) << 1)
                | ((v0 == float(0)) << 2)
                | ((v1 == float(1)) << 3);
            info->failFlag |= failFlag;
        }
#endif      
        return false;
    }

    bool bClip = IsClip(level);
    if (bClip && (IsEps(min, max, eps) || IsLevel(level, max_level))) {
        return testBezierClipL(info, patch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1;
        float tt0 = 0;
        float tt1 = 1;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (GetRangeV(rng, tpatch)) {
                tt0 = rng[0];
                tt1 = rng[1];
                tw = tt1 - tt0;
            }
        }

        if (tw >= 0.4) {
            BezierPatch tmp[2];
            patch.SplitV(tmp, 0.5);
            float vm = (v0 + v1)*0.5;
            float vt[] = { v0,vm,vm,v1 };

            int order[2] = { 0,1 };
            CoarseSort(order, tmp);
            bool bRet = false;
            if (vt[0] > vt[1]) std::swap(vt[0], vt[1]);
            if (vt[2] > vt[3]) std::swap(vt[2], vt[3]);
            for (int i = 0; i < 2; ++i) {
                if (testBezierClipU(info, tmp[order[i]], u0, u1,
                    vt[2 * order[i]], vt[2 * order[i] + 1],
                    zmin, zmax, level + 1, max_level, eps)) {
                    zmax = info.t;
                    bRet = true;
                }
            }
            return bRet;
        }
        else {
            tt0 = std::max<float>(0.0, tt0 - Epsilon::GetUvEps());
            tt1 = std::min<float>(tt1 + Epsilon::GetUvEps(), 1.0);
            float vt[] = { Lerp(v0,v1,tt0),Lerp(v0,v1,tt1) };
            BezierPatch tmp;
            patch.CropV(tmp, tt0, tt1);
            return testBezierClipU(info, tmp, u0, u1,
                vt[0], vt[1], zmin, zmax, level + 1, max_level, eps);
        }
    }
    return false;
}

bool BezierPatchIntersectionImpl::testBezierClipL(UVT& info, const BezierPatch& patch,
    float u0, float u1, float v0, float v1,
    float zmin, float zmax, int level)
{
    float t = 0.0f;
    float u = 0.0f;
    float v = 0.0f;

    if (m_directBilinear) {
        float3 P[4] = {
            patch.Get(  0,   0),
            patch.Get(N-1,   0),
            patch.Get(  0, N-1),
            patch.Get(N-1, N-1),
        };
        if (!m_useTriangle) {
            if (TestBilinearPatch(&t, &u, &v, P, zmin, zmax, m_uvMargin)) {
                u = u0*(1 - u) + u1*u;
                v = v0*(1 - v) + v1*v;
                info.u = u;
                info.v = v;
                info.t = t;
                info.level = level;
                return true;
            }
        }
        else {
            if (TestQuadPlane(&t, &u, &v, P, float3(0.0f), float3(0.0f, 0.0f, 1.0f), zmin, zmax)) {
                u = u0*(1 - u) + u1*u;
                v = v0*(1 - v) + v1*v;
                info.u = u;
                info.v = v;
                info.t = t;
                info.level = level;
                return true;
            }
        }
        int failFlag = ((u0 == float(0)) << 0)
            | ((u1 == float(1)) << 1)
            | ((v0 == float(0)) << 2)
            | ((v1 == float(1)) << 3);

        info.failFlag |= failFlag;
        return false;
    }

    bool bRet = false;
    int nPu = N - 1;
    int nPv = N - 1;

    float du = float(1) / nPu;
    float dv = float(1) / nPv;

    float uu = 0;
    float vv = 0;

    float3 P[4];
    for (int j = 0; j < nPv; ++j) {
        for (int i = 0; i < nPu; ++i) {
            P[0] = patch.Get(  i,   j);
            P[1] = patch.Get(i+1,   j);
            P[2] = patch.Get(  i, j+1);
            P[3] = patch.Get(i+1, j+1);
            if (!m_useTriangle) {
                if (TestBilinearPatch(&t, &u, &v, P, zmin, zmax, m_uvMargin)) {

                    u = Lerp(i*du, (i + 1)*du, u);
                    v = Lerp(j*dv, (j + 1)*dv, v);

                    uu = u;
                    vv = v;

                    zmax = t;

                    u = Lerp(u0, u1, u);
                    v = Lerp(v0, v1, v);

                    info.u = u;
                    info.v = v;
                    info.t = t;
                    info.level = level;
                    bRet = true;
                }
            }
            else {
                if (TestQuadPlane(&t, &u, &v, P, float3(0, 0, 0), float3(0, 0, 1), zmin, zmax)) {
                    u = Lerp(i*du, (i + 1)*du, u);
                    v = Lerp(j*dv, (j + 1)*dv, v);

                    uu = u;
                    vv = v;

                    zmax = t;

                    u = Lerp(u0, u1, u);
                    v = Lerp(v0, v1, v);

                    info.u = u;
                    info.v = v;
                    info.t = t;
                    info.level = level;
                    bRet = true;
                }
            }
        }
    }

    if (bRet) {
        float3 p = patch.Evaluate(float2(uu, vv));
        info.t = p[2];
    }
    else {
        int failFlag = ((u0 == float(0)) << 0)
            | ((u1 == float(1)) << 1)
            | ((v0 == float(0)) << 2)
            | ((v1 == float(1)) << 3);

        info.failFlag |= failFlag;
    }
    return bRet;
}

bool BezierPatchIntersectionImpl::testBezierClipRangeU(UVT& info, const BezierPatch& patch,
    float u0, float u1,
    float v0, float v1, float zmin, float zmax,
    int level, int max_level, float eps)
{
    BezierPatch mpatch(patch, u0, u1, v0, v1);
    BezierPatch tpatch(mpatch);
    RotateU(tpatch);

    float3 min, max;
    tpatch.GetMinMax(min, max, eps*1e-3f);

    if (0.0f < min[0] || max[0] < 0.0f || // x
        0.0f < min[1] || max[1] < 0.0f || // y
        max[2] < zmin || zmax < min[2])   // z
    {
#if USE_MINMAX_FAILFLAG
        int failFlag = ((u0 == float(0)) << 0)
            | ((u1 == float(1)) << 1)
            | ((v0 == float(0)) << 2)
            | ((v1 == float(1)) << 3);
        info->failFlag |= failFlag;
#endif
        return false;
    }

    bool bClip = IsClip(level);
    if (bClip && (IsEps(min, max, eps) || IsLevel(level, max_level))) {
        return testBezierClipL(info, mpatch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1;
        float tt0 = 0;
        float tt1 = 1;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (GetRangeU(rng, tpatch)) {
                tt0 = rng[0];
                tt1 = rng[1];
                tw = tt1 - tt0;
            }
        }

        if (tw >= 0.4) {
            float um = (u0 + u1)*0.5;
            um = float(1.0) - (float(1.0) - um);
            float ut[] = { u0,um,um,u1 };
            int order[2] = { 0,1 };
            bool bRet = false;
            for (int i = 0; i < 2; ++i) {
                if (testBezierClipRangeV(info, patch, ut[2 * order[i]], ut[2 * order[i] + 1],
                    v0, v1, zmin, zmax, level + 1, max_level, eps)) {
                    zmax = info.t;
                    bRet = true;
                }
            }
            return bRet;
        }
        else {
            tt0 = std::max<float>(0.0, tt0 - Epsilon::GetUvEps());
            tt1 = std::min<float>(tt1 + Epsilon::GetUvEps(), 1.0);
            float ut[] = { Lerp(u0,u1,tt0),Lerp(u0,u1,tt1) };
            return testBezierClipRangeV(info, patch, ut[0], ut[1],
                v0, v1, zmin, zmax, level + 1, max_level, eps);
        }
    }
    return false;
}

bool BezierPatchIntersectionImpl::testBezierClipRangeV(UVT& info, const BezierPatch& patch,
    float u0, float u1, float v0, float v1, float zmin, float zmax,
    int level, int max_level, float eps)
{
    BezierPatch mpatch(patch, u0, u1, v0, v1);
    BezierPatch tpatch(mpatch);
    RotateV(tpatch);

    float3 min, max;
    tpatch.GetMinMax(min, max, eps*1e-3f);
    if (0.0f < min[0] || max[0] < 0.0f || // x
        0.0f < min[1] || max[1] < 0.0f || // y
        max[2] < zmin || zmax < min[2])   // z
    {
#if USE_MINMAX_FAILFLAG
        int failFlag = ((u0 == float(0)) << 0)
            | ((u1 == float(1)) << 1)
            | ((v0 == float(0)) << 2)
            | ((v1 == float(1)) << 3);
        info->failFlag |= failFlag;
#endif
        return false;
    }

    bool bClip = IsClip(level);
    if (bClip && (IsEps(min, max, eps) || IsLevel(level, max_level))) {
        return testBezierClipL(info, mpatch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1;
        float tt0 = 0;
        float tt1 = 1;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (GetRangeV(rng, tpatch)) {
                tt0 = rng[0];
                tt1 = rng[1];
                tw = tt1 - tt0;
            }
        }

        if (tw >= 0.4) {
            float vm = (v0 + v1)*0.5;
            vm = float(1.0) - (float(1.0) - vm);
            float vt[] = { v0,vm,vm,v1 };
            int order[2] = { 0,1 };
            bool bRet = false;
            for (int i = 0; i < 2; ++i) {
                if (testBezierClipRangeU(info, patch, u0, u1,
                    vt[2 * order[i]], vt[2 * order[i] + 1],
                    zmin, zmax, level + 1, max_level, eps)) {
                    zmax = info.t;
                    bRet = true;
                }
            }
            return bRet;
        }
        else {
            tt0 = std::max<float>(0.0, tt0 - Epsilon::GetUvEps());
            tt1 = std::min<float>(tt1 + Epsilon::GetUvEps(), 1.0);
            float vt[] = { Lerp(v0,v1,tt0),Lerp(v0,v1,tt1) };
            return testBezierClipRangeU(info, patch, u0, u1,
                vt[0], vt[1], zmin, zmax, level + 1, max_level, eps);
        }
    }
    return false;
}
