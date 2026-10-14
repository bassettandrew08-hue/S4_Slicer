"""
Toolpath-vs-model metrics for tools/model_suite.py (the finished 4-axis G-code against the input STL):

  outside_mm        extruded path more than 1 mm outside the part (plastic in the air: folds, fins)
  surface_rough_pct share of the path near the outer surface (within 0.7 mm) whose segments turn sharply
                    (the path_roughness test, > 30 deg): what shows as a rough wall
  asym_mm           for (nearly) mirror-symmetric parts only: the 90th percentile distance from each extruded
                    point, mirrored across the part's plane, to the nearest extruded point. Bead spacing alone
                    gives ~0.2-0.3 mm; larger means the left and right halves print differently. It sees the union
                    of all layers, so it measures the walls (the part's surface) and Cura's infill pattern (gyroid is
                    not mirror-symmetric: a flat slice of pi 3mm scores 0.67), not the shape of the curved layers
  layer_asym_mm     for the same parts: the 90th percentile of |deformed height at p - deformed height at mirror(p)|
                    over the extruded points p whose mirror is inside the part, i.e. how far (in mm of layer height)
                    a curved layer on one side is from its mirror image on the other. Needs the run's work dir

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


def _centroid(m):
    """Area-weighted surface centroid (exact: random surface samples moved the plane by ~0.02 mm between runs)."""
    V, F = np.asarray(m.vertices), np.asarray(m.triangles)
    A = np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1)
    return (V[F].mean(1) * A[:, None]).sum(0) / A.sum()


def layer_asym(model, stl_path, work_dir, mid):
    """layer_asym_mm (see the module docstring), or None without the work dir's deformed points."""
    import json
    import os
    from s4 import meshio_s4
    dp = os.path.join(work_dir, "deformed_points.npy")
    pu = os.path.join(work_dir, "params_used.json")
    if model not in SYMMETRY or not (os.path.exists(dp) and os.path.exists(pu)):
        return None
    d = json.load(open(pu))["deform"]
    grid = meshio_s4.load_and_tetrahedralize(stl_path, d.get("PART_OFFSET", (0, 0, 0)),
                                             d.get("SURFACE_SIMPLIFY_ERROR", 0.0))
    Dz = np.load(dp)[:, 2]
    if len(Dz) != grid.number_of_points:
        return None
    cells = meshio_s4.cells_of(grid)
    pts = np.asarray(grid.points)
    a = np.radians(SYMMETRY[model])
    n = np.array([np.cos(a), np.sin(a), 0.0])
    c = _centroid(_model(stl_path))
    P = mid[::4]
    Pm = P - 2 * ((P - c) @ n)[:, None] * n

    def height(Q):
        ci = np.asarray(grid.find_containing_cell(Q))
        ok = ci >= 0
        T = pts[cells[np.where(ok, ci, 0)]]
        lam = np.linalg.solve((T[:, 1:] - T[:, :1]).transpose(0, 2, 1), (Q - T[:, 0])[:, :, None])[:, :, 0]
        b = np.c_[1 - lam.sum(1), lam]
        return np.einsum("ij,ij->i", b, Dz[cells[np.where(ok, ci, 0)]]), ok
    h, ok1 = height(P)
    hm, ok2 = height(Pm)
    ok = ok1 & ok2
    return round(float(np.percentile(np.abs(h - hm)[ok], 90)), 3) if ok.any() else None


def measure(model, stl_path, mid, seg_len, rough_deg, work_dir=None):
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
        c = _centroid(m)
        P = mid[::4]
        Pm = P - 2 * ((P - c) @ n)[:, None] * n
        out["asym_mm"] = round(float(np.percentile(cKDTree(mid).query(Pm)[0], 90)), 3)
    else:
        out["asym_mm"] = None
    out["layer_asym_mm"] = layer_asym(model, stl_path, work_dir, mid) if work_dir else None
    return out
