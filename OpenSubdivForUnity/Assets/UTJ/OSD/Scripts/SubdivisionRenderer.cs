using UnityEngine;
using UnityEngine.Rendering;
#if UNITY_EDITOR
using UnityEditor;
#endif


namespace UTJ.OSD
{
    class SubdivisionRenderer : MonoBehaviour
    {
        [SerializeField] Mesh m_mesh;
        [SerializeField] Material[] m_materials;
        [SerializeField] ShadowCastingMode m_shadow = ShadowCastingMode.Off;
        [SerializeField] bool m_receiveShadows = false;
        [SerializeField] LayerSelector m_layer = 0;

        ComputeBuffer m_cbPatches;
#if UNITY_EDITOR
        bool m_dirty = false;
#endif


        public Mesh sharedMesh
        {
            get { return m_mesh; }
            set { m_mesh = value; }
        }

        public Material material
        {
            get { return m_materials != null && m_materials.Length > 0 ? m_materials[0] : null; }
            set { m_materials = new Material[] { value }; }
        }
        public Material[] materials
        {
            get { return m_materials; }
            set { m_materials = value; }
        }



        public void Flush()
        {
            if (m_mesh == null || m_materials == null || m_materials.Length == 0) { return; }

            if (!SystemInfo.supportsComputeShaders) { return; }

            int num_submeshes = System.Math.Min(m_mesh.subMeshCount, m_materials.Length);
            var matrix = GetComponent<Transform>().localToWorldMatrix;
            // issue drawcalls
            for (int si = 0; si < num_submeshes; ++si)
            {
                var material = m_materials[si];
                Graphics.DrawMesh(m_mesh, matrix, material, m_layer, null, si, null, m_shadow, m_receiveShadows);
            }
        }

        public void Release()
        {
            if (m_cbPatches != null) { m_cbPatches.Release(); m_cbPatches = null; }
        }



        void OnDisable()
        {
            Release();
        }

        void LateUpdate()
        {
            Flush();
#if UNITY_EDITOR
            m_dirty = true;
#endif
        }

#if UNITY_EDITOR
        void OnDrawGizmos()
        {
            // force draw particles while paused.
            // using OnDrawGizmos() is dirty workaround but I couldn't find better way...
            if (EditorApplication.isPaused && m_dirty)
            {
                Flush();
                m_dirty = false;
            }
        }
#endif
    }

}