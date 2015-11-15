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
        ComputeBuffer m_buf_bpatches;
        ComputeBuffer m_buf_aabbs;

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
        }
    }
}
