using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Threading;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif


namespace Ist
{
    class BezierPatchContainer : IBezierPatchContainer
    {
        BezierPatchRaw[] m_bpatch_src;
        BezierPatchRaw[] m_bpatch;
        BezierPatchAABB[] m_aabbs;
        bool m_dirty;

        void UpdateBezierPatch()
        {
            if (!m_dirty || m_bpatch_src == null) { return; }
            if (m_bpatch == null)
            {
                m_bpatch = new BezierPatchRaw[m_bpatch_src.Length];
                m_aabbs = new BezierPatchAABB[m_bpatch_src.Length];
            }

            var trans = GetComponent<Transform>().localToWorldMatrix;
            for(int i=0; i<m_bpatch_src.Length; ++i)
            {
                m_bpatch[i] = m_bpatch_src[i];
                m_bpatch[i].Transform(ref trans);
                m_bpatch[i].GetAABB(ref m_aabbs[i]);
            }
        }

        public override BezierPatchRaw[] GetBezierPatches()
        {
            UpdateBezierPatch();
            return m_bpatch;
        }

        public override BezierPatchAABB[] GetAABBs()
        {
            UpdateBezierPatch();
            return m_aabbs;
        }


        public void SetBezierPatches(BezierPatchRaw[] src)
        {
            m_dirty = true;
            m_bpatch_src = src;
        }
    }
}
