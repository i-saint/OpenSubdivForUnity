using UnityEngine;
using System.Collections;
using System.Collections.Generic;


namespace Ist
{
    [AddComponentMenu("MassParticle/GPU Particle/Bezier Patch Collider")]
    [RequireComponent(typeof(IBezierPatchContainer))]
    public class MPGPBezierPatchCollider : MPGPColliderBase
    {
        public Vector3 m_center;
        public Vector3 m_size = Vector3.one;
    
        public override void ActualUpdate()
        {
            var container = GetComponent<IBezierPatchContainer>();
            var bps = container.GetBezierPatches();
            var aabbs = container.GetAABBs();
            MPGPBezierPatchColliderData collider;
            for (int i=0; i<bps.Length; ++i)
            {
                collider.info.aabb.center = aabbs[i].center;
                collider.info.aabb.extents = aabbs[i].extents;
                collider.info.owner_objid = m_id;
                collider.bp = bps[i];
                EachTargets((t) => { t.AddBezierPatchCollider(ref collider); });
            }
        }
    
        void OnDrawGizmos()
        {
            if (!enabled) return;
            Transform t = GetComponent<Transform>();
            Gizmos.color = MPGPImpl.ColliderGizmoColor;
            Gizmos.matrix = t.localToWorldMatrix;
            Gizmos.DrawWireCube(m_center, m_size);
            Gizmos.matrix = Matrix4x4.identity;
        }
    
    }
}