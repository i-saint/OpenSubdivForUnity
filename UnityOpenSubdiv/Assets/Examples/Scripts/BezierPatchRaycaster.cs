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
        public BezierPatchRaw bp;
    };

    public struct TestSplitData
    {
        public BezierPatchRaw bp;
        public BezierPatchRaw dst0;
        public BezierPatchRaw dst1;
        public BezierPatchRaw dst2;
        public BezierPatchRaw dst3;
        public Vector2 uv;
    };

    public struct TestCropData
    {
        public BezierPatchRaw bp;
        public BezierPatchRaw dst;
        public Vector2 uv0;
        public Vector2 uv1;
    };

    public struct TestRaycastData
    {
        public BezierPatchRaw bp;
        public Matrix4x4 bptrans;
        public Vector3 ray_pos;
        public Vector3 ray_dir;
        public float zmin;
        public float zmax;
        public BezierPatchHit hit;
        public Vector3 hit_pos;
        public Vector3 hit_normal;
        public BezierPatchRaw dst;
        public Matrix4x4 zalign;
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
    ComputeBuffer m_buf_split;
    ComputeBuffer m_buf_crop;
    ComputeBuffer m_buf_raycast;
    Vertex[] m_vertices = new Vertex[256];
    TestEvaluateData[] m_data_eval = new TestEvaluateData[1];
    TestSplitData[] m_data_split = new TestSplitData[1];
    TestCropData[] m_data_crop = new TestCropData[1];
    TestRaycastData[] m_data_raycast = new TestRaycastData[1];
    public BezierPatch[] m_splited = new BezierPatch[4];
    public BezierPatch m_croped = new BezierPatch();

    public bool m_use_compute_shader = false;
    public bool m_show_splitted = false;
    public bool m_show_croped = false;


    public void TestEvaluate()
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
            m_bpatch.bpatch.GetRawData(ref m_data_eval[0].bp);

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

    void TestSplit()
    {
        if(m_bpatch == null) { return; }

        if (m_splited[0]==null)
        {
            for(int i=0; i<m_splited.Length; ++i)
            {
                m_splited[i] = new BezierPatch();
            }
        }

        var uv = new Vector2(0.5f, 0.5f);
        if (m_use_compute_shader)
        {
            if (m_buf_split == null)
            {
                m_buf_split = new ComputeBuffer(1, SizeOf<TestSplitData>());
            }
            m_bpatch.bpatch.GetRawData(ref m_data_split[0].bp);
            m_data_split[0].uv = uv;

            int k = m_debug_cs.FindKernel("TestSplit");
            m_buf_split.SetData(m_data_split);
            m_debug_cs.SetBuffer(k, "_TestSplitData", m_buf_split);
            m_debug_cs.Dispatch(k, 1, 1, 1);
            m_buf_split.GetData(m_data_split);

            m_splited[0].SetRawData(ref m_data_split[0].dst0);
            m_splited[1].SetRawData(ref m_data_split[0].dst1);
            m_splited[2].SetRawData(ref m_data_split[0].dst2);
            m_splited[3].SetRawData(ref m_data_split[0].dst3);
        }
        else
        {
            m_bpatch.bpatch.Split(
                ref m_splited[0], ref m_splited[1], ref m_splited[2], ref m_splited[3], ref uv);
        }
    }

    void TestCrop()
    {
        if (m_bpatch == null) { return; }

        var uv0 = new Vector2(0.33f, 0.33f);
        var uv1 = new Vector2(0.66f, 0.66f);

        if (m_use_compute_shader)
        {
            if (m_buf_crop == null)
            {
                m_buf_crop = new ComputeBuffer(1, SizeOf<TestCropData>());
            }
            m_bpatch.bpatch.GetRawData(ref m_data_crop[0].bp);
            m_data_crop[0].uv0 = uv0;
            m_data_crop[0].uv1 = uv1;

            int k = m_debug_cs.FindKernel("TestCrop");
            m_buf_crop.SetData(m_data_crop);
            m_debug_cs.SetBuffer(k, "_TestCropData", m_buf_crop);
            m_debug_cs.Dispatch(k, 1, 1, 1);
            m_buf_crop.GetData(m_data_crop);

            m_croped.SetRawData(ref m_data_crop[0].dst);
        }
        else
        {
            m_bpatch.bpatch.Crop(ref m_croped, ref uv0, ref uv1);
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
        if (m_buf_split != null)
        {
            m_buf_split.Release();
            m_buf_split = null;
        }
        if (m_buf_crop != null)
        {
            m_buf_crop.Release();
            m_buf_crop = null;
        }
        if (m_buf_raycast != null)
        {
            m_buf_raycast.Release();
            m_buf_raycast = null;
        }
    }

    void Update()
    {
    }


    void OnDrawGizmos()
    {
        if (m_bpatch == null) { return; }

        float zmin = 0.0f;
        float zmax = 5.0f;
        float epsilon = 0.01f;
    
        Gizmos.color = Color.red;
        var ray_pos = transform.position;
        var ray_dir = transform.forward;
        var trans = m_bpatch.transform.localToWorldMatrix;
        Gizmos.DrawLine(ray_pos, ray_pos + ray_dir * zmax);

        {
            BezierPatchHit hit = default(BezierPatchHit);
            if (m_bpatch.bpatch.Raycast(ref trans, ray_pos, ray_dir, zmin, zmax, epsilon, ref hit))
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
            m_bpatch.bpatch.GetRawData(ref m_data_raycast[0].bp);
            m_data_raycast[0].bptrans = trans;
            m_data_raycast[0].zmin = 0.0f;
            m_data_raycast[0].zmax = zmax;

            int k = m_debug_cs.FindKernel("TestRaycast");
            m_buf_raycast.SetData(m_data_raycast);
            m_debug_cs.SetBuffer(k, "_TestRaycastData", m_buf_raycast);
            m_debug_cs.Dispatch(k, 1, 1, 1);
            m_buf_raycast.GetData(m_data_raycast);

            Gizmos.color = Color.blue;
            Gizmos.DrawWireSphere(m_data_raycast[0].hit_pos, 0.05f);
        }

        if (m_show_splitted && m_splited[0] != null)
        {
            Color[] colors = new Color[]
            {
                new Color(1.0f, 0.0f, 0.0f),
                new Color(0.0f, 1.0f, 0.0f),
                new Color(0.0f, 0.0f, 1.0f),
                new Color(0.7f, 0.7f, 0.0f),
            };
            for (int i = 0; i < m_splited.Length; ++i)
            {
                m_splited[i].OnDrawGizmo(colors[i]);
            }
        }

        if (m_show_croped)
        {
            m_croped.OnDrawGizmo(new Color(0.7f, 0.0f, 0.7f));
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
        if (m_bpatch == null) { return; }

        if (GUI.Button(new Rect(10, 10, 150, 30), "Split"))
        {
            TestSplit();
        }

        if (GUI.Button(new Rect(10, 50, 150, 30), "Crop"))
        {
            TestCrop();
        }

        if (GUI.Button(new Rect(10, 90, 150, 30), "Evaluate"))
        {
            TestEvaluate();
        }
    }
}
