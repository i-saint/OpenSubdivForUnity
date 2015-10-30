#ifndef BezierPatch_h
#define BezierPatch_h

class BezierPatch
{
public:
    static const int N = 4;
    static const int Ncp = N * N;

    BezierPatch();
    BezierPatch(const float3 *p);
    BezierPatch(const BezierPatch &other);
    BezierPatch(const BezierPatch &other, float u0, float u1, float v0, float v1);
    template<typename MATRIX> BezierPatch(const BezierPatch &other, const MATRIX &m);

    const float3& Get(int i, int j) const;
    void Set(int i, int j, const float3& v);
    float3 Evaluate(float u, float v) const;
    float3 EvaluateDu(float u, float v) const;
    float3 EvaluateDv(float u, float v) const;
    float3 EvaluateNormal(float u, float v) const;

    float3 GetMin() const;
    float3 GetMax() const;
    void GetMinMax(float3 &o_min, float3 &o_max) const;
    void GetMinMax(float3 &o_min, float3 &o_max, float eps) const;

    template<typename MATRIX> BezierPatch& Transform(const MATRIX &m);
    void Split(BezierPatch patches[4], float u, float v) const;
    void SplitU(BezierPatch patches[2], float u) const;
    void SplitV(BezierPatch patches[2], float v) const;
    void CropU(BezierPatch &patch, float u0, float u1) const;
    void CropV(BezierPatch &patch, float v0, float v1) const;
    bool Crop(BezierPatch &patch, float u0, float u1, float v0, float v1) const;
    float3 GetLv() const;
    float3 GetLu() const;
    void Rotate(); // right 90

private:
    void bezierSplit(float3 a[], float3 b[], const float3 p[], float t) const;
    void bezierSplitV(float3 a[], float3 b[], const float3 p[], float t) const;
    void bezierCrop(float3 r[], const float3 p[], float s, float t) const;
    void bezierCropV(float3 r[], const float3 p[], float s, float t) const;
    static void transpose(float3 m[]);
    static float3 evaluate(float t, const float3 * const cp);
    static float3 evaluateD(float t, const float3 * const cp);

private:
    float3 m_cp[Ncp];
};

#endif // !BezierPatch_h
