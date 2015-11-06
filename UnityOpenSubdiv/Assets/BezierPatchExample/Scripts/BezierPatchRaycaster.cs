using System.Collections;
using UnityEngine;
using Ist;

[ExecuteInEditMode]
public class BezierPatchRaycaster : MonoBehaviour
{
    public BezierPatchEditor m_bpatch;

    void OnDrawGizmos()
    {
        float max_distance = 10.0f;

        Gizmos.color = Color.red;
        var pos = transform.position;
        var dir = transform.forward;
        Gizmos.DrawLine(pos, pos + dir * max_distance);

        if(m_bpatch != null)
        {
            BezierPatchHit hit = default(BezierPatchHit);
            if(m_bpatch.bpatch.Raycast(pos, dir, max_distance, ref hit))
            {
                //var hit_pos = pos + dir * hit.t;
                var hit_pos = m_bpatch.bpatch.Evaluate(hit.uv);
                Gizmos.DrawWireCube(hit_pos, Vector3.one*0.1f);
            }
        }
    }
}
