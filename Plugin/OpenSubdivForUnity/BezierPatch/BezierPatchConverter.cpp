#include "pch.h"
#include "osuInternal.h"
#include "BezierPatch.h"
#include "BezierPatchConverter.h"

using namespace OpenSubdiv::Osd;

BezierPatchConverter::BezierPatchConverter(int num_vertices)
    : m_vertex_buffer(nullptr)
    , m_topology_refiner(nullptr)
    , m_patch_table(nullptr)
    , m_stencil_table(nullptr)
{
    m_vertex_buffer = CpuVertexBuffer::Create(3, num_vertices);

    bool use_single_crease_patch = false;
    int max_iteration_level = 1;

    // topology refiner
    {
        Shape * shape = nullptr;
        // todo: build shape

        OpenSubdiv::Sdc::SchemeType type = GetSdcType(*shape);
        OpenSubdiv::Sdc::Options options = GetSdcOptions(*shape);
        OpenSubdiv::Far::TopologyRefinerFactory<Shape>::Options opt(type, options);
        m_topology_refiner = OpenSubdiv::Far::TopologyRefinerFactory<Shape>::Create(*shape, opt);
    }
    // refine
    {
        OpenSubdiv::Far::TopologyRefiner::AdaptiveOptions options(max_iteration_level);
        options.useSingleCreasePatch = use_single_crease_patch;
        m_topology_refiner->RefineAdaptive(options);
    }

    // patch table
    {
        OpenSubdiv::Far::PatchTableFactory::Options options;
        options.useSingleCreasePatch = use_single_crease_patch;
        options.endCapType = OpenSubdiv::Far::PatchTableFactory::Options::ENDCAP_BSPLINE_BASIS;
        options.maxIsolationLevel = max_iteration_level;
        m_patch_table = OpenSubdiv::Far::PatchTableFactory::Create(*m_topology_refiner, options);

    }

    // stencil table
    {
        OpenSubdiv::Far::StencilTableFactory::Options options;
        options.generateOffsets = true;
        options.generateIntermediateLevels = true;
        m_stencil_table = OpenSubdiv::Far::StencilTableFactory::Create(*m_topology_refiner, options);
    }
}

BezierPatchConverter::~BezierPatchConverter()
{
    delete m_vertex_buffer;
    m_vertex_buffer = nullptr;

    delete m_topology_refiner;
    m_topology_refiner = nullptr;
}

void BezierPatchConverter::Convert(int num_vertices, const float3 *vertices)
{
    m_vertex_buffer->UpdateData(&vertices[0][0], 0, num_vertices);

    int num_control_vertices = m_stencil_table->GetNumControlVertices();
    OpenSubdiv::Osd::BufferDescriptor src_desc(0, 3, 3);
    OpenSubdiv::Osd::BufferDescriptor dst_desc(3 * num_control_vertices, 3, 3);
    OpenSubdiv::Osd::CpuEvaluator::EvalStencils(
        m_vertex_buffer, src_desc,
        m_vertex_buffer, dst_desc,
        m_stencil_table);
}

const std::vector<BezierPatch>& BezierPatchConverter::GetBezierPatches() const
{
    return m_bezier_patches;
}
