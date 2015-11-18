using System.Collections;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif
using Ist;


[ExecuteInEditMode]
public class BezierPatchRaycaster : MonoBehaviour
{
    public struct Vertex
    {
        public Vector3 vertex;
        public Vector3 normal;
    }

    public struct TestEvaluateData
    {
        public BezierPatchRaw bpatch;
    };

    public struct TestRaycastData
    {
        public Vector3 ray_pos;
        public Vector3 ray_dir;
        public BezierPatchRaw bpatch;
        public Matrix4x4 bptrans;
        public float zmin;
        public float zmax;
        public BezierPatchHit hit;
        public Vector3 hit_pos;
        public Vector3 hit_normal;
    };

    public static int SizeOf<T>()
    {
        return System.Runtime.InteropServices.Marshal.SizeOf(typeof(T));
    }


    public BezierPatchEditor m_bpatch;
    public ComputeShader m_debug_cs;
    public Mesh m_mesh;
    ComputeBuffer m_buf_vertices;
    ComputeBuffer m_buf_eval;
    ComputeBuffer m_buf_raycast;
    Vertex[] m_vertices = new Vertex[256];
    TestEvaluateData[] m_data_eval = new TestEvaluateData[1];
    TestRaycastData[] m_data_raycast = new TestRaycastData[1];
    public BezierPatch[] m_splited = new BezierPatch[4];

    public void UpdatePreviewMesh()
    {
        const int div = 16;
        const int divsq = div * div;
        bool update_indices = false;

        if (m_mesh == null)
        {
            update_indices = true;

            GameObject go = new GameObject();
            go.name = "Bezier Patch Mesh";
            go.GetComponent<Transform>().SetParent(GetComponent<Transform>());

            var mesh_filter = go.AddComponent<MeshFilter>();
            var mesh_renderer = go.AddComponent<MeshRenderer>();
#if UNITY_EDITOR
            mesh_renderer.sharedMaterial = AssetDatabase.LoadAssetAtPath<Material>("Assets/Examples/Materials/Default.mat");
#endif

            m_mesh = new Mesh();
            mesh_filter.sharedMesh = m_mesh;
        }

        // update vertices
        if (m_debug_cs != null)
        {
            if (m_buf_eval == null)
            {
                m_buf_eval = new ComputeBuffer(1, 12*16);
                m_buf_vertices = new ComputeBuffer(256, 24);
            }
            m_data_eval[0].bpatch = m_bpatch.GetBezierPatches()[0];

            int k = m_debug_cs.FindKernel("TestEvaluate");
            m_buf_eval.SetData(m_data_eval);
            m_debug_cs.SetBuffer(k, "_Vertices", m_buf_vertices);
            m_debug_cs.SetBuffer(k, "_TestEvaluateData", m_buf_eval);
            m_debug_cs.Dispatch(k, divsq, 1, 1);
            m_buf_vertices.GetData(m_vertices);
            m_buf_eval.GetData(m_data_eval);

            {
                var vertices = new Vector3[divsq];
                var normals = new Vector3[divsq];
                for (int i = 0; i < divsq; ++i)
                {
                    vertices[i] = m_vertices[i].vertex;
                    normals[i] = m_vertices[i].normal;
                }
                m_mesh.vertices = vertices;
                m_mesh.normals = normals;
            }


            if (update_indices)
            {
                var indices = new int[divsq * 6];
                for (int y = 0; y < div - 1; ++y)
                {
                    for (int x = 0; x < div - 1; ++x)
                    {
                        indices[(y * div + x) * 6 + 0] = (y + 0) * div + (x + 0);
                        indices[(y * div + x) * 6 + 1] = (y + 1) * div + (x + 0);
                        indices[(y * div + x) * 6 + 2] = (y + 1) * div + (x + 1);

                        indices[(y * div + x) * 6 + 3] = (y + 0) * div + (x + 0);
                        indices[(y * div + x) * 6 + 4] = (y + 1) * div + (x + 1);
                        indices[(y * div + x) * 6 + 5] = (y + 0) * div + (x + 1);
                    }
                }
                m_mesh.SetIndices(indices, MeshTopology.Triangles, 0);
            }
        }
    }

    void Split()
    {
        if (m_splited[0]==null)
        {
            for(int i=0; i<m_splited.Length; ++i)
            {
                m_splited[i] = new BezierPatch();
            }
        }
        if(m_bpatch != null)
        {
            var uv = new Vector2(0.5f, 0.5f);
            m_bpatch.bpatch.Split(
                ref m_splited[0], ref m_splited[1], ref m_splited[2], ref m_splited[3], ref uv);
        }
    }



    void OnDisable()
    {
        if (m_buf_eval != null)
        {
            m_buf_eval.Release();
            m_buf_eval = null;
            m_buf_vertices.Release();
            m_buf_vertices = null;
        }
        if (m_buf_raycast != null)
        {
            m_buf_raycast.Release();
            m_buf_raycast = null;
        }
    }

    void Update()
    {
        UpdatePreviewMesh();
    }


    void OnDrawGizmos()
    {
        if (m_bpatch == null) { return; }

        float max_distance = 5.0f;
    
        Gizmos.color = Color.red;
        var ray_pos = transform.position;
        var ray_dir = transform.forward;
        var trans = m_bpatch.transform.localToWorldMatrix;
        Gizmos.DrawLine(ray_pos, ray_pos + ray_dir * max_distance);

        {
            BezierPatchHit hit = default(BezierPatchHit);
            if (m_bpatch.bpatch.Raycast(ref trans, ray_pos, ray_dir, max_distance, ref hit))
            {
                var hit_pos = ray_pos + ray_dir * hit.t;
                //var hit_pos = m_bpatch.bpatch.Evaluate(hit.uv); // same as above
                Gizmos.DrawWireSphere(hit_pos, 0.05f);
            }
        }

        if(m_debug_cs != null)
        {
            if(m_buf_raycast == null)
            {
                m_buf_raycast = new ComputeBuffer(1, SizeOf<TestRaycastData>());
            }
            m_data_raycast[0].ray_pos = ray_pos;
            m_data_raycast[0].ray_dir = ray_dir;
            m_data_raycast[0].bpatch = m_bpatch.GetBezierPatches()[0];
            m_data_raycast[0].bptrans = trans;
            m_data_raycast[0].zmin = 0.0f;
            m_data_raycast[0].zmax = max_distance;

            int k = m_debug_cs.FindKernel("TestRaycast");
            m_buf_raycast.SetData(m_data_raycast);
            m_debug_cs.SetBuffer(k, "_TestRaycastData", m_buf_raycast);
            m_debug_cs.Dispatch(k, 1, 1, 1);
            m_buf_raycast.GetData(m_data_raycast);

            Gizmos.color = Color.blue;
            Gizmos.DrawWireSphere(m_data_raycast[0].hit_pos, 0.05f);
        }

        if (m_splited[0] != null)
        {
            Color[] colors = new Color[]
            {
                new Color(1.0f, 0.0f, 0.0f),
                new Color(0.0f, 1.0f, 0.0f),
                new Color(0.0f, 0.0f, 0.0f),
                new Color(0.6f, 0.6f, 0.0f),
            };
            for (int i = 0; i < m_splited.Length; ++i)
            {
                m_splited[i].OnDrawGizmo(colors[i]);
            }
        }
    }

    // for debug
    /*
    void Raycast()
    {
        if (m_bpatch != null)
        {
            var pos = transform.position;
            var dir = transform.forward;
            float max_distance = 10.0f;
            BezierPatchHit hit = default(BezierPatchHit);
            bool ret = m_bpatch.bpatch.Raycast(pos, dir, max_distance, ref hit);
        }
    }
        */

    void OnGUI()
    {
        if (GUI.Button(new Rect(10, 10, 150, 30), "Split"))
        {
            Split();
        }
    }
}
