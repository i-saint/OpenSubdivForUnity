#include "pch.h"
#include "BezierPatch.h"
#include "BezierPatchConverter.h"

using namespace OpenSubdiv::Osd;

BezierPatchConverter::BezierPatchConverter(int num_vertices)
    : m_vertex_buffer(nullptr)
    , m_topology_refiner(nullptr)
    , m_stencil_table(nullptr)
{
    m_vertex_buffer = CpuVertexBuffer::Create(3, num_vertices);

    {
        Shape * shape = nullptr;
        // todo: build shape

        OpenSubdiv::Sdc::SchemeType type = GetSdcType(*shape);
        OpenSubdiv::Sdc::Options options = GetSdcOptions(*shape);
        OpenSubdiv::Far::TopologyRefinerFactory<Shape>::Options opt(type, options);
        m_topology_refiner = OpenSubdiv::Far::TopologyRefinerFactory<Shape>::Create(*shape, opt);
    }
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
