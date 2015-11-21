#ifndef BezierPatchConverter_h
#define BezierPatchConverter_h

class BezierPatchConverter
{
public:
    BezierPatchConverter(int num_vertices);
    ~BezierPatchConverter();
    void Clear();
    void Convert(int num_vertices, const float3 *vertices);
    void Convert(int num_indices, const float3 *vertices, const uint32_t *indices);
    void Convert(int num_indices, const float3 *vertices, const uint16_t *indices);

    const std::vector<BezierPatch>& GetBezierPatches() const;

private:
    std::vector<BezierPatch> m_bezier_patches;
    std::vector<float> m_sharpnesses;

    OpenSubdiv::Osd::CpuVertexBuffer *m_vertex_buffer;
    OpenSubdiv::Far::TopologyRefiner *m_topology_refiner;
    OpenSubdiv::Far::PatchParamTable m_patch_param_table;
    const OpenSubdiv::Far::StencilTable *m_stencil_table;
};

#endif // BezierPatchConverter_h
