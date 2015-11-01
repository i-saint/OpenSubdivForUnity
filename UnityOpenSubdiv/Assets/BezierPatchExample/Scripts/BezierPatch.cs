using UnityEngine;
using System.Collections;

[ExecuteInEditMode]
public class BezierPatch : MonoBehaviour
{
    Vector3[] m_cp = new Vector3[16];
    Transform[] m_cpobj;
    ComputeBuffer m_cb;

    void OnEnable()
    {
        if(m_cpobj == null)
        {
            var trans = GetComponent<Transform>();
            m_cpobj = new Transform[16];
            float span = 2.0f / 3.0f;
            for (int y = 0; y < 4; ++y)
            {
                for (int x = 0; x < 4; ++x)
                {
                    var go = new GameObject();
                    go.name = "Control Point ["+ y +"][" + x +"]";
                    go.AddComponent<BezierPatchControlPoint>();
                    var t = go.GetComponent<Transform>();
                    t.position = new Vector3(-1.0f + span * x, 0.0f, -1.0f + span * y);
                    t.SetParent(trans);
                    m_cpobj[y * 4 + x] = t;
                }
            }

        }
    }

    void OnDisable()
    {

    }

    void LateUpdate()
    {
        for(int i=0; i< m_cp.Length; ++i)
        {
            m_cp[i] = m_cpobj[i].position;
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
