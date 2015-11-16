using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Threading;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif


namespace Ist
{
    [RequireComponent(typeof(MeshRenderer))]
    [RequireComponent(typeof(IBezierPatchContainer))]
    class BezierPatchRenderer : MonoBehaviour
    {
        [SerializeField] Mesh m_bound_mesh;
        [SerializeField] Material m_material;
        Material m_material_copy; // I need to copy material because MaterialPropertyBlock don't have SetBuffer()
        ComputeBuffer m_buf_bpatches;
        ComputeBuffer m_buf_aabbs;


#if UNITY_EDITOR
        void Reset()
        {
        }
#endif

        void OnDestroy()
        {
            if (m_buf_bpatches != null)
            {
                m_buf_bpatches.Release();
                m_buf_bpatches = null;
            }
            if (m_buf_aabbs != null)
            {
                m_buf_aabbs.Release();
                m_buf_aabbs = null;
            }
        }

        void OnPreRender()
        {
            if(m_material_copy == null)
            {
                m_material_copy = new Material(m_material);
            }

            var cont = GetComponent<IBezierPatchContainer>();
            var bpatches = cont.GetBezierPatches();
            var aabbs = cont.GetAABBs();
            if (m_buf_bpatches == null)
            {
                m_buf_bpatches = new ComputeBuffer(bpatches.Length, BezierPatchRaw.size);
                m_material_copy.SetBuffer("_BezierPatches", m_buf_bpatches);
            }
            if (m_buf_aabbs == null)
            {
                m_buf_aabbs = new ComputeBuffer(aabbs.Length, BezierPatchAABB.size);
                m_material_copy.SetBuffer("_AABBs", m_buf_aabbs);
            }
            m_buf_bpatches.SetData(bpatches);
            m_buf_aabbs.SetData(aabbs);
        }
    }
}
