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
        public const int size = 12 * 16;

        public Vector3
            cp00, cp01, cp02, cp03,
            cp04, cp05, cp06, cp07,
            cp08, cp09, cp10, cp11,
            cp12, cp13, cp14, cp15;
    }

    [System.Serializable]
    public struct BezierPatchAABB
    {
        public const int size = 12 * 2;

        public Vector3 center, extents;
    }

    [System.Serializable]
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


        public void GetRawData(ref BezierPatchRaw dst)
        {
            // I HATE C#
            dst.cp00 = cp[ 0]; dst.cp01 = cp[ 1]; dst.cp02 = cp[ 2]; dst.cp03 = cp[ 3];
            dst.cp04 = cp[ 4]; dst.cp05 = cp[ 5]; dst.cp06 = cp[ 6]; dst.cp07 = cp[ 7];
            dst.cp08 = cp[ 8]; dst.cp09 = cp[ 9]; dst.cp10 = cp[10]; dst.cp11 = cp[11];
            dst.cp12 = cp[12]; dst.cp13 = cp[13]; dst.cp14 = cp[11]; dst.cp15 = cp[15];
        }

        public void GetAABB(ref BezierPatchAABB dst)
        {
            Vector3 min = cp[0];
            Vector3 max = cp[0];
            for(int i=1; i<cp.Length; ++i)
            {
                min = Vector3.Min(min, cp[i]);
                max = Vector3.Max(max, cp[i]);
            }
            dst.center = (max + min) * 0.5f;
            dst.extents = (max - min) * 0.5f;
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
