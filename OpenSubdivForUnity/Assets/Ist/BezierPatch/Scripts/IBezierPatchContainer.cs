
using UnityEngine;

namespace Ist
{
    public abstract class IBezierPatchContainer : MonoBehaviour
    {
        public abstract BezierPatchRaw[] GetBezierPatches();
        public abstract BezierPatchAABB[] GetAABBs();
    }
}
