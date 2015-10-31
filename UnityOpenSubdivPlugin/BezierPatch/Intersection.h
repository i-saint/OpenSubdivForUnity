#ifndef Intersection_h
#define Intersection_h

#include "BezierPatch.h"


struct Ray
{
    float3 org;
    float3 dir;
    float3 invDir;
    int dirSign[3];
    int depth;
    bool hasDifferential;

    float3 dDdx;
    float3 dDdy;
};

struct Curve
{
    static const int LENGTH = 4;
    float2 operator[](int i) const { return v[i]; }
    float2 &operator[](int i) { return v[i]; }
    float2 v[LENGTH];
};

struct Intersection
{
    float t;
    float u;
    float v;
    uint32_t faceID;
    uint32_t matID;

    // patch info
    uint32_t patchID;
    uint32_t level;
    uint32_t clipLevel;

    // for measurements
    float eps;
    uint32_t maxLevel;

    // triangle
    uint32_t f0;
    uint32_t f1;
    uint32_t f2;

    float3 position;
    float3 geometricNormal;
    float3 normal;
    float3 tangent;
    float3 binormal;
    float2 texcoord;
};

class BezierPatchIntersection
{
public:
    static const int N = 4;
    static const int  DEFAULT_MAX_LEVEL = 10;
    static const int  MAX_ITERATION = 1000;

    struct RangeAABB {
        float tmin, tmax;
    };
    struct UVT {
        UVT() : u(0), v(0), t(0), level(0) { }
        float u, v, t;
        int level;
        int failFlag;
    };


    BezierPatchIntersection(const BezierPatch& patch);
    ~BezierPatchIntersection();

    void SetEpsilon(float eps);
    void SetMaxLevel(int maxLevel);
    void SetUVMergin(float uvMargin);
    void SetCropUV(bool cropUV);
    void SetUseBezierClip(bool useBezierClip);
    void SetUseTriangle(bool useTriangle);
    void SetDirectBilinear(bool directBilinear);
    void SetWatertightFlag(int wcpFlag);

    bool Test(Intersection* info, const Ray& r, float tmin, float tmax);
    float ComputeEpsilon(Ray const & r, float eps) const;

protected:
    bool testInternal(Intersection* info, const Ray& r, float tmin, float tmax);
    bool testBezierPatch(UVT* info, BezierPatch const & patch, float zmin, float zmax, float eps);


    bool testBezierClipU(UVT* info, BezierPatch const & patch,
        float u0, float u1, float v0, float v1,
        float zmin, float zmax,
        int level, int max_level, float eps);

    bool testBezierClipV(UVT *info, const BezierPatch& patch,
        float u0, float u1, float v0, float v1, float zmin, float zmax,
        int level, int max_level, float eps);

    bool testBezierClipL(UVT* info, const BezierPatch& patch,
        float u0, float u1, float v0, float v1,
        float zmin, float zmax, int level);

    bool testBezierClipRangeU(UVT* info, const BezierPatch& patch,
        float u0, float u1,
        float v0, float v1, float zmin, float zmax,
        int level, int max_level, float eps);

    bool testBezierClipRangeV(UVT *info, const BezierPatch& patch,
        float u0, float u1, float v0, float v1, float zmin, float zmax,
        int level, int max_level, float eps);

    template<typename Ray>
    bool intersectAABB(RangeAABB *rng, float3 const & min, float3 const & max,
        Ray const & r, float tmin, float tmax) const;

    BezierPatch m_patch;
    float m_uRange[2];
    float m_vRange[2];
    float3 m_min;
    float3 m_max;
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

#endif // Intersection_h
