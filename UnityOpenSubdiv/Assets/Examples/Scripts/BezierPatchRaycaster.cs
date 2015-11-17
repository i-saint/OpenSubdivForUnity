using System.Collections;
using UnityEngine;
using Ist;


[ExecuteInEditMode]
public class BezierPatchRaycaster : MonoBehaviour
{
    public struct DebugData
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
    ComputeBuffer m_cb;
    DebugData[] m_data = new DebugData[1];

    void OnDisable()
    {
        if(m_cb != null)
        {
            m_cb.Release();
            m_cb = null;
        }
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
            if(m_cb == null)
            {
                m_cb = new ComputeBuffer(1, SizeOf<DebugData>());
            }
            m_data[0].ray_pos = ray_pos;
            m_data[0].ray_dir = ray_dir;
            m_data[0].bpatch = m_bpatch.GetBezierPatches()[0];
            m_data[0].bptrans = trans;
            m_data[0].zmin = 0.0f;
            m_data[0].zmax = max_distance;
            m_cb.SetData(m_data);
            m_debug_cs.Dispatch(0, 1, 0, 0);
            m_cb.GetData(m_data);

            Gizmos.color = Color.blue;
            Gizmos.DrawWireSphere(m_data[0].hit_pos, 0.05f);
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
    
    void OnGUI()
    {
        if (GUI.Button(new Rect(10, 10, 150, 30), "Raycast"))
        {
            Raycast();
        }
    }
    */
}
