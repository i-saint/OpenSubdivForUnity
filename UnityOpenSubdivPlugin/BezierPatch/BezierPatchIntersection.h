#ifndef BezierPatchIntersection_h
#define BezierPatchIntersection_h

#include "BezierPatch.h"


struct Ray
{
    float3 org;
    float3 dir;
};

struct BezierPatchHit
{
    float t;
    float u;
    float v;
    uint32_t clip_level;
};

bool BezierPatchRayIntersection(const BezierPatch &bp, const Ray &ray, BezierPatchHit &hit);


class BezierPatchIntersectionImpl
{
public:
    static const int N = 4;
    static const int  DEFAULT_MAX_LEVEL = 10;
    static const int  MAX_ITERATION = 1000;


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

    bool Test(BezierPatchHit &info, const Ray& r, float tmin, float tmax);

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

#endif // BezierPatchIntersection_h
