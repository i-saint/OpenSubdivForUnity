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
    BezierPatch(const BezierPatch &other, const float3x3 &m);
    BezierPatch(const BezierPatch &other, const float4x4 &m);

    float3& operator[](int i) { return m_cp[i]; }
    const float3& operator[](int i) const { return m_cp[i]; }
    const float3& Get(int i, int j) const;
    void Set(int i, int j, const float3& v);

    float3 GetMin() const;
    float3 GetMax() const;
    void GetMinMax(float3 &o_min, float3 &o_max, float epsilon=0.0f) const;

    float3 Evaluate(const float2& uv) const;
    float3 EvaluateDu(const float2& uv) const;
    float3 EvaluateDv(const float2& uv) const;
    float3 EvaluateNormal(const float2& uv) const;
    void Split(BezierPatch dst[4], float u, float v) const;
    void SplitU(BezierPatch &dst0, BezierPatch &dst1, float u) const;
    void SplitV(BezierPatch &dst0, BezierPatch &dst1, float v) const;
    void CropU(BezierPatch &dst, float u0, float u1) const;
    void CropV(BezierPatch &dst, float v0, float v1) const;
    void Crop(BezierPatch &dst, float u0, float u1, float v0, float v1) const;
    float3 GetLv() const;
    float3 GetLu() const;

    void Transpose();
    void Transform(const float3x3 &m);
    void Transform(const float4x4 &m);

    // F: void =>[](float3 &control_point);
    template<class F>
    void ForEach(const F &f)
    {
        for (int i = 0; i < Ncp; ++i) { f(m_cp[i]); }
    }

private:
    float3 m_cp[Ncp];
};





#endif // BezierPatch_h
