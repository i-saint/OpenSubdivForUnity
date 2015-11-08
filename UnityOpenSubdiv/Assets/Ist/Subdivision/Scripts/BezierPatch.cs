using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif

namespace Ist
{
    // struct for store ComputeBuffer
    [System.Serializable]
    public struct BezierPatchRaw
    {
        Vector3
            cp00, cp01, cp02, cp03,
            cp10, cp11, cp12, cp13,
            cp20, cp21, cp22, cp23,
            cp30, cp31, cp32, cp33;
    }

    public struct BezierPatchHit
    {
        public Vector2 uv;
        public float t;
        public uint clip_level;
    };



    [System.Serializable]
    public class BezierPatch
    {
        public Vector3[] cp = DefaultControlPoints();


        public Vector3 Evaluate(Vector2 uv)
        {
            return uosBezierPatchEvaluate(ref cp[0], ref uv);
        }

        public Vector3 EvaluateNormal(Vector2 uv)
        {
            return uosBezierPatchEvaluateNormal(ref cp[0], ref uv);
        }

        public bool Raycast(Vector3 orig, Vector3 dir, float max_distance, ref BezierPatchHit hit)
        {
            return uosBezierPatchRaycast(ref cp[0], ref orig, ref dir, max_distance, ref hit);
        }

        public bool Raycast(ref Matrix4x4 trans, Vector3 orig, Vector3 dir, float max_distance, ref BezierPatchHit hit)
        {
            return uosBezierPatchRaycastWithTransform(ref cp[0], ref trans, ref orig, ref dir, max_distance, ref hit);
        }


        #region impl
        static Vector3[] DefaultControlPoints()
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


        [DllImport("UnityOpenSubdiv")]
        static extern Vector3 uosBezierPatchEvaluate(ref Vector3 bpatch, ref Vector2 uv);

        [DllImport("UnityOpenSubdiv")]
        static extern Vector3 uosBezierPatchEvaluateNormal(ref Vector3 bpatch, ref Vector2 uv);

        [DllImport("UnityOpenSubdiv")]
        static extern bool uosBezierPatchRaycast(ref Vector3 bpatch, ref Vector3 orig, ref Vector3 dir, float max_distance, ref BezierPatchHit hit);

        [DllImport("UnityOpenSubdiv")]
        static extern bool uosBezierPatchRaycastWithTransform(ref Vector3 bpatch, ref Matrix4x4 bptrans, ref Vector3 orig, ref Vector3 dir, float max_distance, ref BezierPatchHit hit);
        #endregion
    }
}
