using UnityEngine;
using System.Collections;

[ExecuteInEditMode]
public class BezierPatchEditor : MonoBehaviour
{
    public bool m_lock;
    public Vector3[] m_cp = InitializeBezierPatch();
    Transform[] m_cpobj;
    ComputeBuffer m_cb;
    

    static Vector3[] InitializeBezierPatch()
    {
        var cp = new Vector3[16];
        float span = 2.0f / 3.0f;
        for (int y = 0; y < 4; ++y)
        {
            for (int x = 0; x < 4; ++x)
            {
                cp[4 * y + x] = new Vector3(-1.0f + span * x, 0.0f, -1.0f + span * y);
            }
        }
        return cp;
    }

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
                    t.position = m_cp[i];
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
            for (int i = 0; i < m_cp.Length; ++i)
            {
                m_cp[i] = m_cpobj[i].position;
            }
        }
    }


    void OnDrawGizmos()
    {
        Gizmos.color = Color.cyan;
        for (int y = 0; y < 4; ++y)
        {
            for (int x = 0; x < 3; ++x)
            {
                Gizmos.DrawLine(m_cp[y*4 + x], m_cp[y*4 + x+1]);
            }
        }
        for (int y = 0; y < 3; ++y)
        {
            for (int x = 0; x < 4; ++x)
            {
                Gizmos.DrawLine(m_cp[y * 4 + x], m_cp[(y+1) * 4 + x]);
            }
        }
    }
}
