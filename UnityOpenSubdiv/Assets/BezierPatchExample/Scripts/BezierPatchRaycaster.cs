using System.Collections;
using UnityEngine;
using Ist;

[ExecuteInEditMode]
public class BezierPatchRaycaster : MonoBehaviour
{
    public BezierPatchEditor m_bpatch;

    void OnDrawGizmos()
    {
        Gizmos.color = Color.red;
        var pos = transform.position;
        var dir = transform.forward;
        Gizmos.DrawLine(pos, pos + dir*10.0f);

        if(m_bpatch != null)
        {
            BezierPatchHit hit = default(BezierPatchHit);
            float max_distance = 10.0f;
            if(m_bpatch.bpatch.Raycast(pos, dir, max_distance, ref hit))
            {
                Gizmos.DrawWireCube(pos + dir * hit.t, Vector3.one*0.1f);
                //Gizmos.DrawWireCube(m_bpatch.bpatch.Evaluate(hit.uv), Vector3.one * 0.1f);
            }
        }
    }
}
