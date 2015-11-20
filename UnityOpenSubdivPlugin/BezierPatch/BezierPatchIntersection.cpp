#include "pch.h"
#include "BezierPatch.h"
#include "BezierPatchIntersection.h"

#define trace(...)
#define EPSILON 1e-4f
#define UVEPSILON (1.0f / 32.0f)
#define DEFAULT_MAX_LEVEL 10
#define MAX_ITERATION 1000



class BezierPatchIntersectionImpl
{
public:
    BezierPatchIntersectionImpl(const BezierPatch& patch);
    ~BezierPatchIntersectionImpl();

    void SetEpsilon(float eps);
    void SetMaxLevel(int maxLevel);
    void SetUVMergin(float uvMargin);
    void SetCropUV(bool cropUV);
    void SetUseBezierClip(bool useBezierClip);
    void SetUseTriangle(bool useTriangle);
    void SetDirectBilinear(bool directBilinear);
    void SetWatertightFlag(int wcpFlag);

    bool Test(BezierPatchHit& info, const Ray& r, float tmin, float tmax);

private:
    struct UVT
    {
        UVT() : u(0), v(0), t(0), level(0) { }
        float u, v, t;
        int level;
        int failFlag;
    };

    bool testInternal(BezierPatchHit& info, const Ray& r, float tmin, float tmax);
    bool testBezierPatch(UVT& info, BezierPatch const & patch, float zmin, float zmax, float eps);


    bool testBezierClipU(UVT& info, BezierPatch const & patch,
        float u0, float u1, float v0, float v1,
        float zmin, float zmax,
        int level, int max_level, float eps);

    bool testBezierClipV(UVT& info, const BezierPatch& patch,
        float u0, float u1, float v0, float v1, float zmin, float zmax,
        int level, int max_level, float eps);

    bool testBezierClipL(UVT& info, const BezierPatch& patch,
        float u0, float u1, float v0, float v1,
        float zmin, float zmax, int level);

    bool testBezierClipRangeU(UVT& info, const BezierPatch& patch,
        float u0, float u1,
        float v0, float v1, float zmin, float zmax,
        int level, int max_level, float eps);

    bool testBezierClipRangeV(UVT& info, const BezierPatch& patch,
        float u0, float u1, float v0, float v1, float zmin, float zmax,
        int level, int max_level, float eps);


    const BezierPatch &m_patch;
    float m_eps;
    int m_maxLevel;
    float m_uvMargin;
    bool m_cropUV;
    bool m_useBezierClip;
    bool m_useTriangle;
    bool m_directBilinear;
    int m_wcpFlag;

    int m_count;
};


static inline float3x3 BPIRotate2D_(const float3 &dx)
{
    float2 x = normalize(float2(dx.x, dx.y));
    float2 y = float2(-x[1], x[0]);
    return float3x3(
        x[0], x[1], 0.0f,
        y[0], y[1], 0.0f,
        0.0f, 0.0f, 1.0f
    );
}


static inline float4x4 BPIZAlign_(const float3& p, const float3& dir)
{
    float3 z = dir;
    int plane = 0;
    if (fabs(z[1]) < fabs(z[plane])) plane = 1;
    if (fabs(z[2]) < fabs(z[plane])) plane = 2;

    float3 up = (plane == 0) ? float3(1.0f, 0.0f, 0.0f) :
               ((plane == 1) ? float3(0.0f, 1.0f, 0.0f) :
                               float3(0.0f, 0.0f, 1.0f));
    float3 y = normalize(cross(z, up));
    float3 x = cross(y, z);

    float4x4 rot = float4x4(
        x[0], y[0], z[0], 0.0f,
        x[1], y[1], z[1], 0.0f,
        x[2], y[2], z[2], 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f);
    float4x4 trs = float4x4(
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
       -p[0],-p[1],-p[2], 1.0f);
    return rot * trs;
}



static inline void BPIRotateU_(BezierPatch& patch)
{
    float3 dx = patch.GetLv();
    patch.Transform(BPIRotate2D_(dx));
}

static inline void BPIRotateV_(BezierPatch& patch)
{
    float3 dx = patch.GetLu();
    patch.Transform(BPIRotate2D_(dx));
}

static inline bool BPIIsEps_(const float3& min, const float3& max, float eps)
{
    //return ((max[1]-min[1])<=eps);
    float xw = max[0] - min[0];
    float yw = max[1] - min[1];
    return (xw <= eps && yw <= eps);
}

static inline bool BPIIsLevel_(int level, int max_level)
{
    return (level >= max_level);
}

static inline bool BPIIsClip_(int level)
{
    static const int div_level[] =
    {
        1,
        1,1,
        2,2,//4
        3,3,3,3,//8
        4,4,4,4,4,4,4,4,//16
    };
    int nlevel = div_level[4];
    return nlevel * 2 <= level;
}

// --------------------------------------------------------------------

struct BPCurve
{
    const float2& operator[](int i) const { return v[i]; }
    float2& operator[](int i) { return v[i]; }

    float2 v[4];
};

static inline bool BPIScanMinMax_(const BPCurve& p)
{
    int nUpper = 0, nLower = 0;
    for (int i = 0; i < 4; ++i) {
        if (p[i][1] > 0) nUpper++;
        else nLower++;
        if (nUpper && nLower) return true;
    }
    return false;
}

static inline float BPISlope_(const float2& a, const float2& b)
{
    float d0 = b.x - a.x;
    float d1 = b.y - a.y;
    return fabs(d0 / d1);
}

static inline float BPIDist_(const float2 &p0, const float2& p1)
{
    float a = fabs(p0.y);
    float b = fabs(p1.y);
    float h = a + b;
    float t = p0.x + a*(p1.x - p0.x) / h;
    return t;
}

static inline int BPIScanConvex_(float ts[], const BPCurve & p)
{
    if (!BPIScanMinMax_(p)) {
        return 0;
    }
    int n = 0;
    {
        int k = -1;
        int l = -1;
        float current = std::numeric_limits<float>::max();
        for (int i = 1; i < 4; ++i) {
            if (p[i][1] * p[0][1] < 0) {
                float s = BPISlope_(p[i], p[0]);
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
                    float s = BPISlope_(p[i], p[k]);
                    if (current < s) {
                        current = s;
                        l = i;
                    }
                }
            }
            if (l >= 0) {
                ts[n++] = BPIDist_(p[l], p[k]);
            }
        }
    }
    {
        int k = -1;
        int l = -1;
        float current = std::numeric_limits<float>::max();
        for (int i = 0; i < 4 - 1; ++i) {
            if (p[i][1] * p[4 - 1][1] < 0) {
                float s = BPISlope_(p[i], p[4 - 1]);
                if (s < current) {
                    current = s;
                    k = i;
                }
            }
        }
        if (k >= 0) {
            current = 0;
            for (int i = k + 1; i < 4; ++i) {
                if (p[i][1] * p[k][1] < 0) {
                    float s = BPISlope_(p[i], p[k]);
                    if (current < s) {
                        current = s;
                        l = i;
                    }
                }
            }
            if (l >= 0) {
                ts[n++] = BPIDist_(p[k], p[l]);
            }
        }
    }
    return n;
}

static inline bool BPIXCheck_(float rng[2], const BPCurve& curve)
{
    float t[4] = { 0 };
    int nn = BPIScanConvex_(t, curve);
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

static inline bool BPIGetRangeU_(float out[2], const BezierPatch& bp)
{
    float t0 = 1.0f;
    float t1 = 0.0f;
    BPCurve curve;
    float delta = 1.0 / (4 - 1);
    for (int i = 0; i < 4; ++i) {
        curve[i][0] = i*delta;
    }
    float rng[2];
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            curve[i][1] = (bp.Get(i, j))[1];
        }
        if (BPIXCheck_(rng, curve)) {
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

static inline bool BPIGetRangeV_(float out[2], const BezierPatch& bp)
{
    float t0 = 1;
    float t1 = 0;
    BPCurve curve;
    float delta = 1.0 / (4 - 1);
    for (int i = 0; i < 4; ++i) {
        curve[i][0] = i*delta;
    }
    float rng[2];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            curve[j][1] = (bp.Get(i, j))[1];
        }
        if (BPIXCheck_(rng, curve)) {
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


static inline bool BPISolveBilinearPatch_(float* t, float* u, float* v, float tmin, float tmax, float tt, float uu, float vv)
{
    if (tmin <= tt && tt <= tmax) {
        *t = tt;
        *u = uu;
        *v = vv;
        return true;
    }
    return false;
}

static inline int BPISolve2_(float root[2], float const coeff[3])
{
    float A = coeff[0];
    float B = coeff[1];
    float C = coeff[2];
    if (fabs(A) <= std::numeric_limits<float>::min()) {
        if (fabs(B) <= std::numeric_limits<float>::min()) return 0;
        float x = -C / B;
        root[0] = x;
        return 1;
    }
    else {
        float D = B*B - 4 * A*C;
        if (D < 0) {
            return 0;
        }
        else if (D == 0) {
            float x = -0.5*B / A;
            root[0] = x;
            return 1;
        }
        else {
            float x1 = (fabs(B) + sqrt(D)) / (2.0 * A);
            if (B >= 0) {
                x1 = -x1;
            }
            float x2 = C / (A * x1);
            if (x1 > x2) std::swap(x1, x2);
            root[0] = x1;
            root[1] = x2;
            return 2;
        }
    }
}

static inline float BPIComputeU_(float A1, float A2, float B1, float B2, float C1, float C2, float D1, float D2, float v)
{
    //return div((v*(C1-C2)+(D1-D2)),(v*(A2-A1)+(B2-B1)));
    float a = v*A2 + B2;
    float b = v*(A2 - A1) + B2 - B1;
    if (fabs(b) >= fabs(a)) {
        return (v*(C1 - C2) + D1 - D2) / b;
    }
    else {
        return (-v*C2 - D2) / a;
    }
}

static inline float BPIComputeT_(float a, float b, float c, float d, float iq, float u, float v)
{
    return ((u*v)*a + u*b + v*c + d)*iq;
}

static inline bool BPITestBilinearPatch_(float *t, float *u, float *v, const float3 p[4], float tmin, float tmax, float uvMargin)
{
    trace("Test bilinear (%f, %f, %f) (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n",
        p[0][0], p[0][1], p[0][2],
        p[1][0], p[1][1], p[1][2],
        p[2][0], p[2][1], p[2][2],
        p[3][0], p[3][1], p[3][2]);

    const float3& p00 = p[0];
    const float3& p10 = p[1];
    const float3& p01 = p[2];
    const float3& p11 = p[3];

    static const int nPlane = 2;
    float3 a = p11 - p10 - p01 + p00;
    float3 b = p10 - p00;
    float3 c = p01 - p00;
    float3 d = p00;

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
    int nRet = BPISolve2_(root, coeff);

    if (nRet) {
        bool bRet = false;
        trace("Solve bilinear:");
        for (int i = 0; i < nRet; ++i) {
            float vv = root[i];
            trace(" vv = %f, ", vv);
            if (0 - uvMargin <= vv && vv <= 1 + uvMargin) {//TODO
                vv = std::max(float(0), std::min(vv, float(1)));
                float uu = BPIComputeU_(A1, A2, B1, B2, C1, C2, D1, D2, vv);
                trace(" uu = %f, ", uu);
                if (0 - uvMargin <= uu && uu <= 1 + uvMargin) {//TODO
                    uu = std::max(float(0), std::min(uu, float(1)));
                    float tt = BPIComputeT_(a[nPlane], b[nPlane], c[nPlane], d[nPlane],
                        float(1), uu, vv);
                    trace("uv=(%f, %f)", uu, vv);
                    if (BPISolveBilinearPatch_(t, u, v, tmin, tmax, tt, uu, vv)) {
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

static inline bool BPITestTriangle_(float *tt, float *uu, float *vv,
    const float3& p0, const float3& p1, const float3& p2, const float3& org, const float3& dir, float tmin, float tmax)
{
    float t, u, v;
    float3 e1(p0 - p1);
    float3 e2(p0 - p2);
    float3 bDir(cross(e2, dir));
    float iM = dot(e1, bDir);

    if (0 < iM) {
        //p0-org
        float3 vOrg(p0 - org);

        u = dot(vOrg, bDir);
        if (u < 0 || iM < u) return false;

        float3 vE(cross(e1, vOrg));

        v = dot(dir, vE);
        if (v < 0 || iM < u + v) return false;

        t = -dot(e2, vE);
    }
    else if (iM < 0) {
        //p0-org
        float3 vOrg(p0 - org);//JKL

        u = dot(vOrg, bDir);
        if (u > 0 || iM > u) return false;

        float3 vE(cross(e1, vOrg));

        v = dot(dir, vE);
        if (v > 0 || iM > u + v) return false;

        t = -dot(e2, vE);
    }
    else {
        return false;
    }

    iM = 1.0f / iM;
    t *= iM;
    if (t < tmin || tmax < t) return false;
    u *= iM;
    v *= iM;

    *tt = t;
    *uu = u;
    *vv = v;

    return true;
}


static inline bool BPITestQuadPlane_(float *t, float *u, float *v,
    const float3 p[4], const float3& org, const float3& dir, float tmin, float tmax)
{
    const float3& p00 = p[0];
    const float3& p10 = p[1];
    const float3& p01 = p[2];
    const float3& p11 = p[3];

    float tt, uu, vv;
    bool bRet = false;
    if (BPITestTriangle_(&tt, &uu, &vv, p00, p01, p10, org, dir, tmin, tmax))
    {
        float ww = float(1) - (uu + vv);
        *u = ww*float(0) + uu*float(0) + vv*float(1);//00 - 01 - 10
        *v = ww*float(0) + uu*float(1) + vv*float(0);//00 - 01 - 10
        *t = tt;
        tmax = tt;
        bRet = true;
    }
    if (BPITestTriangle_(&tt, &uu, &vv, p10, p01, p11, org, dir, tmin, tmax))
    {
        float ww = float(1) - (uu + vv);
        *u = ww*float(1) + uu*float(0) + vv*float(1);//10 - 01 - 11
        *v = ww*float(0) + uu*float(1) + vv*float(1);//10 - 01 - 11
        *t = tt;
        tmax = tt;
        bRet = true;
    }

    return bRet;
}


BezierPatchIntersectionImpl::BezierPatchIntersectionImpl(const BezierPatch &patch) :
    m_patch(patch),
    m_eps(EPSILON),
    m_maxLevel(DEFAULT_MAX_LEVEL),
    m_uvMargin(true),
    m_cropUV(true), m_useBezierClip(true),
    m_useTriangle(false),
    m_directBilinear(false),
    m_wcpFlag(0)
{
    m_count = 0;
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


bool BezierPatchIntersectionImpl::Test(BezierPatchHit& info, const Ray& r, float tmin, float tmax)
{
    return testInternal(info, r, tmin, tmax);
}


bool BezierPatchIntersectionImpl::testInternal(BezierPatchHit& info, const Ray& r, float tmin, float tmax)
{
    float4x4 mat = BPIZAlign_(r.orig, r.dir); // getZAlign
    BezierPatch patch(m_patch, mat);

    bool bRet = false;
    UVT uvt;
    uvt.failFlag = 0;

    if (testBezierPatch(uvt, patch, tmin, tmax, m_eps)) {
        float t = uvt.t;
        float u = uvt.u;
        float v = uvt.v;

        info.t = t;
        info.uv.x = u;
        info.uv.y = v;
        info.clip_level = uvt.level;
        tmax = t;
        bRet = true;
        trace("hit t = %f, uv = (%f, %f)\n", t, u, v);
    }

    if (uvt.failFlag == 0 || m_wcpFlag == 0) {
        return bRet;
    }

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
            0.0f, 0.0f,
            0.0f, 0.5f,
            0.5f, 0.0f,
            0.5f, 0.5f
        };
        BezierPatch tmp[2];
        BezierPatch children[4];
        m_patch.SplitU(tmp[0], tmp[1], 0.5f);
        tmp[0].SplitV(children[0], children[1], 0.5f);
        tmp[1].SplitV(children[2], children[3], 0.5f);

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

            info.t = t;
            info.uv.x = u;
            info.uv.y = v;
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
    if (0.0f < min.x || max.x < 0.0f || //x
        0.0f < min.y || max.y < 0.0f || //y
        max.z < zmin || zmax < min.z )  //z
    {
        return false;
    }

    if (m_cropUV) return testBezierClipRangeU(info, patch, 0.0f, 1.0f, 0.0f, 1.0f, zmin, zmax, 0, m_maxLevel, eps);
    else          return testBezierClipU(info, patch, 0.0f, 1.0f, 0.0f, 1.0f, zmin, zmax, 0, m_maxLevel, eps);
}



bool BezierPatchIntersectionImpl::testBezierClipU(UVT& info, const BezierPatch & patch,
    float u0, float u1, float v0, float v1,
    float zmin, float zmax,
    int level, int max_level, float eps)
{
    BezierPatch tpatch(patch);
    BPIRotateU_(tpatch);

    if (++m_count > MAX_ITERATION || u0 > u1 || v0 > v1) {
        return false;
    }

    float3 min, max;
    patch.GetMinMax(min, max, 1e-3f);

    if (0.0f < min.x || max.x < 0.0f || // x
        0.0f < min.y || max.y < 0.0f || // y
        max.z < zmin || zmax < min.z )  // z
    {
        return false;
    }

    bool bClip = BPIIsClip_(level);
    if (bClip && (BPIIsEps_(min, max, eps) || BPIIsLevel_(level, max_level))) {
        return testBezierClipL(info, patch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1;
        float tt0 = 0;
        float tt1 = 1;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (BPIGetRangeU_(rng, tpatch)) {
                tt0 = rng[0];
                tt1 = rng[1];
                tw = tt1 - tt0;
            }
        }

        if (tw >= 0.4) {
            BezierPatch tmp[2];
            patch.SplitU(tmp[0], tmp[1], 0.5f);
            float um = (u0 + u1)*0.5f;
            float ut[] = { u0,um,um,u1 };

            int order[2] = { 0,1 };
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
            tt0 = std::max<float>(0.0, tt0 - UVEPSILON);
            tt1 = std::min<float>(tt1 + UVEPSILON, 1.0);
            float ut[] = { glm::mix(u0,u1,tt0), glm::mix(u0,u1,tt1) };
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
    BPIRotateV_(tpatch);

    if (++m_count > MAX_ITERATION || u0 > u1 || v0 > v1) {
        return false;
    }

    trace("testBezierClipV (%f, %f) - (%f, %f) z:%f, %f  level=%d\n",
        u0, u1, v0, v1, zmin, zmax, level);

    float3 min, max;
    tpatch.GetMinMax(min, max, 1e-3f);

    if (0.0f < min[0] || max[0] < 0.0f || // x
        0.0f < min[1] || max[1] < 0.0f || // y
        max[2] < zmin || zmax < min[2]) // z
    {
        return false;
    }

    bool bClip = BPIIsClip_(level);
    if (bClip && (BPIIsEps_(min, max, eps) || BPIIsLevel_(level, max_level))) {
        return testBezierClipL(info, patch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1;
        float tt0 = 0;
        float tt1 = 1;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (BPIGetRangeV_(rng, tpatch)) {
                tt0 = rng[0];
                tt1 = rng[1];
                tw = tt1 - tt0;
            }
        }

        if (tw >= 0.4) {
            BezierPatch tmp[2];
            patch.SplitV(tmp[0], tmp[1], 0.5f);
            float vm = (v0 + v1)*0.5f;
            float vt[] = { v0,vm,vm,v1 };

            int order[2] = { 0,1 };
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
            tt0 = std::max<float>(0.0f, tt0 - UVEPSILON);
            tt1 = std::min<float>(tt1 + UVEPSILON, 1.0f);
            float vt[] = { glm::mix(v0,v1,tt0), glm::mix(v0,v1,tt1) };
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
            patch.Get(4-1,   0),
            patch.Get(  0, 4-1),
            patch.Get(4-1, 4-1),
        };
        if (!m_useTriangle) {
            if (BPITestBilinearPatch_(&t, &u, &v, P, zmin, zmax, m_uvMargin)) {
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
            if (BPITestQuadPlane_(&t, &u, &v, P, float3(0.0f), float3(0.0f, 0.0f, 1.0f), zmin, zmax)) {
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
    int nPu = 4 - 1;
    int nPv = 4 - 1;

    float du = 1.0f / nPu;
    float dv = 1.0f / nPv;

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
                if (BPITestBilinearPatch_(&t, &u, &v, P, zmin, zmax, m_uvMargin)) {

                    u = glm::mix(i*du, (i + 1)*du, u);
                    v = glm::mix(j*dv, (j + 1)*dv, v);

                    uu = u;
                    vv = v;

                    zmax = t;

                    u = glm::mix(u0, u1, u);
                    v = glm::mix(v0, v1, v);

                    info.u = u;
                    info.v = v;
                    info.t = t;
                    info.level = level;
                    bRet = true;
                }
            }
            else {
                if (BPITestQuadPlane_(&t, &u, &v, P, float3(0, 0, 0), float3(0, 0, 1), zmin, zmax)) {
                    u = glm::mix(i*du, (i + 1)*du, u);
                    v = glm::mix(j*dv, (j + 1)*dv, v);

                    uu = u;
                    vv = v;

                    zmax = t;

                    u = glm::mix(u0, u1, u);
                    v = glm::mix(v0, v1, v);

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
    float u0, float u1, float v0, float v1, float zmin, float zmax,
    int level, int max_level, float eps)
{
    BezierPatch mpatch(patch, u0, u1, v0, v1);
    BezierPatch tpatch(mpatch);
    BPIRotateU_(tpatch);

    float3 min, max;
    tpatch.GetMinMax(min, max, eps*1e-3f);

    if (0.0f < min[0] || max[0] < 0.0f || // x
        0.0f < min[1] || max[1] < 0.0f || // y
        max[2] < zmin || zmax < min[2])   // z
    {
        return false;
    }

    bool bClip = BPIIsClip_(level);
    if (bClip && (BPIIsEps_(min, max, eps) || BPIIsLevel_(level, max_level))) {
        return testBezierClipL(info, mpatch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1.0f;
        float tt0 = 0.0f;
        float tt1 = 1.0f;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (BPIGetRangeU_(rng, tpatch)) {
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
            tt0 = std::max<float>(0.0, tt0 - UVEPSILON);
            tt1 = std::min<float>(tt1 + UVEPSILON, 1.0);
            float ut[] = { glm::mix(u0,u1,tt0), glm::mix(u0,u1,tt1) };
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
    BPIRotateV_(tpatch);

    float3 min, max;
    tpatch.GetMinMax(min, max, eps*1e-3f);
    if (0.0f < min[0] || max[0] < 0.0f || // x
        0.0f < min[1] || max[1] < 0.0f || // y
        max[2] < zmin || zmax < min[2])   // z
    {
        return false;
    }

    bool bClip = BPIIsClip_(level);
    if (bClip && (BPIIsEps_(min, max, eps) || BPIIsLevel_(level, max_level))) {
        return testBezierClipL(info, mpatch, u0, u1, v0, v1, zmin, zmax, level);
    }
    else {
        float tw = 1;
        float tt0 = 0;
        float tt1 = 1;

        if (m_useBezierClip & bClip) {
            float rng[2];
            if (BPIGetRangeV_(rng, tpatch)) {
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
            tt0 = std::max<float>(0.0, tt0 - UVEPSILON);
            tt1 = std::min<float>(tt1 + UVEPSILON, 1.0);
            float vt[] = { glm::mix(v0,v1,tt0), glm::mix(v0,v1,tt1) };
            return testBezierClipRangeU(info, patch, u0, u1,
                vt[0], vt[1], zmin, zmax, level + 1, max_level, eps);
        }
    }
    return false;
}




bool BezierPatchRaycast(const BezierPatch &bp, const Ray &ray, float zmax, BezierPatchHit &hit)
{
    BezierPatchIntersectionImpl impl(bp);
    return impl.Test(hit, ray, 0.0f, zmax);
}

bool BezierPatchRaycast(const BezierPatch &bp, const float4x4 &bptrans, const Ray &ray, float zmax, BezierPatchHit &hit)
{
    BezierPatchIntersectionImpl impl(BezierPatch(bp, bptrans));
    return impl.Test(hit, ray, 0.0f, zmax);
}


// debugging compute shader 
/*
#ifndef BPI_MAX_STACK_DEPTH
#define BPI_MAX_STACK_DEPTH 20
#endif
#ifndef BPI_MAX_LOOP
#define BPI_MAX_LOOP 1000
#endif
#ifndef BPI_EPS
#define BPI_EPS 0.1
#endif

struct BPIWorkingBuffer
{
    BezierPatch source; // input
    BezierPatch crop;
    BezierPatch rotate;
    float4 uv_range; // input
};

bool BPITriangleIntersect_(
    float& tout, float& uout, float& vout,
    float3 p0, float3 p1, float3 p2,
    float3 ray_org, float3 ray_dir)
{
    float3 e1, e2;
    float3 p, s, q;

    e1 = p1 - p0;
    e2 = p2 - p0;
    p = cross(ray_dir, e2);

    float det = dot(e1, p);
    float inv_det = 1.0 / det;

    s = ray_org - p0;
    q = cross(s, e1);

    float u = dot(s, p) * inv_det;
    float v = dot(q, ray_dir) * inv_det;
    float t = dot(e2, q) * inv_det;

    if (u < 0.0 || u > 1.0) return false;
    if (v < 0.0 || u + v > 1.0) return false;
    if (t < 0.0 || t > tout) return false;

    tout = t;
    uout = u;
    vout = v;
    return true;
}

bool BPITestBezierClipL_(const BezierPatch &bp, float3 &uvt, const float2 &uv0, const float2 &uv1, float zmin, float zmax)
{
    // TODO (NO_DIRECT)
    // DIRECT_BILINEAR
    float3 p0, p1, p2, p3;
    float3 ray_org = float3(0.0, 0.0, 0.0);
    float3 ray_dir = float3(0.0, 0.0, 1.0);
    p0 = bp[0];
    p1 = bp[3];
    p2 = bp[12];
    p3 = bp[15];
    bool ret = false;
    float t = zmax, uu = 0.0, vv = 0.0;
    if (BPITriangleIntersect_(t, uu, vv, p0, p2, p1, ray_org, ray_dir)) {
        float ww = 1.0 - (uu + vv);
        float u = ww*0.0 + uu*0.0 + vv*1.0; //00 - 01 - 10
        float v = ww*0.0 + uu*1.0 + vv*0.0; //00 - 01 - 10
        uvt.x = glm::mix(uv0.x, uv1.x, u);
        uvt.y = glm::mix(uv0.y, uv1.y, v);
        uvt.z = t;
        ret = true;
    }
    if (BPITriangleIntersect_(t, uu, vv, p1, p2, p3, ray_org, ray_dir)) {
        float ww = 1.0 - (uu + vv);
        float u = ww*1.0 + uu*0.0 + vv*1.0; //10 - 01 - 11
        float v = ww*0.0 + uu*1.0 + vv*1.0; //10 - 01 - 11
        uvt.x = glm::mix(uv0.x, uv1.x, u);
        uvt.y = glm::mix(uv0.y, uv1.y, v);
        uvt.z = t;
        ret = true;
    }
    return ret;
}

bool BPITestBezierPatch_(BPIWorkingBuffer &work, BezierPatchHit &info, float zmin, float zmax, float eps)
{
    // non-recursive iteration
    float4 range_stack[BPI_MAX_STACK_DEPTH];
    int stack_index = 0;
    range_stack[0] = work.uv_range;

    int i;
    for (i = 0; i < BPI_MAX_LOOP && stack_index >= 0; ++i) {

        // pop a patch range and crop
        float u0 = range_stack[stack_index].x;
        float u1 = range_stack[stack_index].y;
        float v0 = range_stack[stack_index].z;
        float v1 = range_stack[stack_index].w;
        --stack_index;

        work.source.Crop(work.crop, float2(u0, v0), float2(u1, v1));
        float3 LU = work.crop[3] - work.crop[0];
        float3 LV = work.crop[12] - work.crop[0];
        bool clipU = length(LU) > length(LV);

        float3 bmin, bmax;
        // rotate and bmin/bmax
        float3 dx = clipU
            ? work.crop[12] - work.crop[0] + work.crop[15] - work.crop[3]
            : work.crop[3] - work.crop[0] + work.crop[15] - work.crop[12];
        work.rotate = work.crop;
        work.rotate.Transform(BPIRotate2D_(dx));
        work.rotate.GetMinMax(bmin, bmax, eps*1e-3);

        // out
        if (0.0 < bmin.x || bmax.x < 0.0 || 0.0 < bmin.y || bmax.y < 0.0 || bmax.z < zmin || zmax < bmin.z) {
            continue;
        }

        // if it's small enough, test bilinear.
        if ((bmax.x - bmin.x) < eps || (bmax.y - bmin.y) < eps) {
            float3 uvt;
            if (BPITestBezierClipL_(work.crop, uvt, float2(u0, v0), float2(u1, v1), zmin, zmax)) {
                info.uv.x = uvt.x;
                info.uv.y = uvt.y;
                info.t = uvt.z;
                zmax = info.t;
                return true;
            }
            // find another intersection
            continue;
        }
        info.clip_level = i;

        // push children ranges
        if (clipU) {
            float um = (u0 + u1)*0.5;
            range_stack[++stack_index] = float4(u0, um, v0, v1);
            range_stack[++stack_index] = float4(um, u1, v0, v1);
        }
        else {
            float vm = (v0 + v1)*0.5;
            range_stack[++stack_index] = float4(u0, u1, v0, vm);
            range_stack[++stack_index] = float4(u0, u1, vm, v1);
        }

        if (stack_index >= BPI_MAX_STACK_DEPTH - 1) break;
    }
    return false;
}

bool BezierPatchRaycast(const BezierPatch &bp, const Ray &ray, float zmax, BezierPatchHit &hit)
{
    BPIWorkingBuffer work;
    work.source = bp;
    work.source.Transform(BPIZAlign_(ray.orig, ray.dir));
    work.uv_range = float4(0.0f, 1.0f, 0.0f, 1.0f);
    return BPITestBezierPatch_(work, hit, 0.0f, zmax, BPI_EPS);
}

bool BezierPatchRaycast(const BezierPatch &bp, const float4x4 &bptrans, const Ray &ray, float zmax, BezierPatchHit &hit)
{
    BPIWorkingBuffer work;
    work.source = bp;
    work.source.Transform(bptrans);
    work.source.Transform(BPIZAlign_(ray.orig, ray.dir));
    work.uv_range = float4(0.0f, 1.0f, 0.0f, 1.0f);
    return BPITestBezierPatch_(work, hit, 0.0f, zmax, BPI_EPS);
}
*/
