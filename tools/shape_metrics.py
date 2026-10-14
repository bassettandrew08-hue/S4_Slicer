"""
Toolpath-vs-model metrics for tools/model_suite.py (the finished 4-axis G-code against the input STL):

  outside_mm        extruded path more than 1 mm outside the part (plastic in the air: folds, fins)
  surface_rough_pct share of the path near the outer surface (within 0.7 mm) whose segments turn sharply
                    (the path_roughness test, > 30 deg): what shows as a rough wall
  asym_mm           for (nearly) mirror-symmetric parts only: the 90th percentile distance from each extruded
                    point, mirrored across the part's plane, to the nearest extruded point. Bead spacing alone
                    gives ~0.2-0.3 mm; larger means the left and right halves print differently

The toolpath is in the bed frame, where the part is centred on its bounding box in x and y with z = 0 at its
bottom (meshio_s4.load_and_tetrahedralize), so the STL is moved the same way here.
"""
import numpy as np
import open3d as o3d
from scipy.spatial import cKDTree

# mirror planes found by tools: (normal angle in the xy plane, deg) through the part's surface centroid;
# benchy and dino are only nearly symmetric (text, small details)
SYMMETRY = {"pi 3mm": 90.0, "benchy upsidedown tilted": 90.0, "dino": 3.0}


def _model(stl_path):
    m = o3d.io.read_triangle_mesh(stl_path)
    V = np.asarray(m.vertices)
    shift = np.array([(V[:, 0].min() + V[:, 0].max()) / 2, (V[:, 1].min() + V[:, 1].max()) / 2, V[:, 2].min()])
    m.vertices = o3d.utility.Vector3dVector(V - shift)
    return m


def measure(model, stl_path, mid, seg_len, rough_deg):
    """mid, seg_len, rough_deg: per bead segment (quality.path_roughness's mid and rough, plus lengths)."""
    m = _model(stl_path)
    scene = o3d.t.geometry.RaycastingScene()
    scene.add_triangles(o3d.t.geometry.TriangleMesh.from_legacy(m))
    sd = scene.compute_signed_distance(o3d.core.Tensor(mid.astype(np.float32))).numpy()
    out = {"outside_mm": round(float(seg_len[sd > 1.0].sum()), 1)}
    near = np.abs(sd) < 0.7
    out["surface_rough_pct"] = round(100.0 * float((rough_deg[near] > 30).mean()), 2) if near.any() else None
    if model in SYMMETRY:
        a = np.radians(SYMMETRY[model])
        n = np.array([np.cos(a), np.sin(a), 0.0])
        c = np.asarray(m.sample_points_uniformly(100000).points).mean(0)
        P = mid[::4]
        Pm = P - 2 * ((P - c) @ n)[:, None] * n
        out["asym_mm"] = round(float(np.percentile(cKDTree(mid).query(Pm)[0], 90)), 3)
    else:
        out["asym_mm"] = None
    return out
