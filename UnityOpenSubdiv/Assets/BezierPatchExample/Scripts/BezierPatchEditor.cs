using System.Collections;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif
using Ist;

[ExecuteInEditMode]
public class BezierPatchEditor : MonoBehaviour
{
    public BezierPatch m_bpatch = new BezierPatch();
    public bool m_lock;
    Transform[] m_cpobj;
    ComputeBuffer m_cb;


#if UNITY_EDITOR
    public void GenerateMesh()
    {
        GameObject go = new GameObject();
        go.name = "Bezier Patch Mesh";
        var mesh_filter = go.AddComponent<MeshFilter>();
        var mesh_renderer = go.AddComponent<MeshRenderer>();
        mesh_renderer.sharedMaterial = AssetDatabase.LoadAssetAtPath<Material>("Assets/BezierPatchExample/Materials/Default.mat");

        var mesh = new Mesh();
        {
            const int div = 16;
            const int divsq = div * div;
            var indices = new int[divsq * 6];
            var vertices = new Vector3[divsq];
            var normals = new Vector3[divsq];
            var span = new Vector2(1.0f / (div-1), 1.0f / (div-1));

            for (int y = 0; y < div; ++y)
            {
                for (int x = 0; x < div; ++x)
                {
                    int i = y * div + x;
                    var uv = new Vector2(span.x*x, span.y*y);
                    vertices[i] = m_bpatch.Evaluate(uv);
                    normals[i] = -m_bpatch.EvaluateNormal(uv);
                }
            }
            for (int y = 0; y < div-1; ++y)
            {
                for (int x = 0; x < div-1; ++x)
                {
                    indices[(y * div + x) * 6 + 0] = (y + 0) * div + (x + 0);
                    indices[(y * div + x) * 6 + 1] = (y + 1) * div + (x + 0);
                    indices[(y * div + x) * 6 + 2] = (y + 1) * div + (x + 1);

                    indices[(y * div + x) * 6 + 3] = (y + 0) * div + (x + 0);
                    indices[(y * div + x) * 6 + 4] = (y + 1) * div + (x + 1);
                    indices[(y * div + x) * 6 + 5] = (y + 0) * div + (x + 1);
                }
            }
            mesh.vertices = vertices;
            mesh.normals = normals;
            mesh.SetIndices(indices, MeshTopology.Triangles, 0);
        }

        mesh_filter.sharedMesh = mesh;
    }
#endif

    void DestroyControlPointes()
    {
        if (m_cpobj != null)
        {
            for (int i = 0; i < m_cpobj.Length; ++i)
            {
                if(m_cpobj[i] != null)
                {
                    DestroyImmediate(m_cpobj[i].gameObject);
                }
            }
            m_cpobj = null;
        }
    }

    void ConstructControlPoints()
    {
        if (m_cpobj == null || m_cpobj.Length != 16)
        {
            m_cpobj = new Transform[16];
        }

        var trans = GetComponent<Transform>();
        for (int y = 0; y < 4; ++y)
        {
            for (int x = 0; x < 4; ++x)
            {
                int i = y * 4 + x;
                if (m_cpobj[i] == null)
                {
                    var go = new GameObject();
                    go.name = "Control Point [" + y + "][" + x + "]";
                    go.AddComponent<BezierPatchControlPoint>();
                    var t = go.GetComponent<Transform>();
                    t.position = m_bpatch.cp[i];
                    t.SetParent(trans);
                    m_cpobj[i] = t;
                }
            }
        }
    }


    void OnDisable()
    {
        DestroyControlPointes();
    }

    void Update()
    {
        if(m_lock)
        {
            DestroyControlPointes();
        }
        else
        {
            ConstructControlPoints();
            for (int i = 0; i < m_bpatch.cp.Length; ++i)
            {
                m_bpatch.cp[i] = m_cpobj[i].position;
            }
        }
    }


    void OnDrawGizmos()
    {
        Gizmos.color = Color.cyan;
        var cp = m_bpatch.cp;
        for (int y = 0; y < 4; ++y)
        {
            for (int x = 0; x < 3; ++x)
            {
                Gizmos.DrawLine(cp[y*4 + x], cp[y*4 + x+1]);
            }
        }
        for (int y = 0; y < 3; ++y)
        {
            for (int x = 0; x < 4; ++x)
            {
                Gizmos.DrawLine(cp[y * 4 + x], cp[(y+1) * 4 + x]);
            }
        }
    }
}
