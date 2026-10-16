"""
Toolpath-vs-model metrics for tools/model_suite.py (the finished 4-axis G-code against the input STL):

  outside_mm        extruded path more than 1 mm outside the part (plastic in the air: folds, fins)
  surface_rough_pct share of the path near the outer surface (within 0.7 mm) whose segments turn sharply
                    (the path_roughness test, > 30 deg): what shows as a rough wall
  asym_mm           for (nearly) mirror-symmetric parts only: the 90th percentile distance from each extruded
                    point, mirrored across the part's plane, to the nearest extruded point. Bead spacing alone
                    gives ~0.2-0.3 mm; larger means the left and right halves print differently

The toolpath is in the bed frame, where the part is centred on its bounding box in x and y with z = 0 at its
bottom, then moved by -PART_OFFSET (meshio_s4.load_and_tetrahedralize), so the STL is moved the same way here
(part_offset_of reads the run's PART_OFFSET).
"""
import json
import os
import re

import numpy as np
import open3d as o3d
from scipy.spatial import cKDTree

# mirror planes found by tools: (normal angle in the xy plane, deg) through the part's surface centroid;
# benchy and dino are only nearly symmetric (text, small details)
SYMMETRY = {"pi 3mm": 90.0, "benchy upsidedown tilted": 90.0, "dino": 3.0}


def part_offset_of(gcode_path=None, work_dir=None):
    """The run's PART_OFFSET: from the G-code's settings header ("; s4: PART_OFFSET = x, y, z"), else from the work
    dir's params_used.json, else 0."""
    if gcode_path and os.path.exists(gcode_path):
        with open(gcode_path, errors="replace") as fh:
            for line in fh:
                m = re.match(r"; s4: PART_OFFSET = (.*)", line)
                if m:
                    return np.array([float(x) for x in re.split(r"[,\s\[\]]+", m.group(1).strip()) if x], float)
                if line.startswith("; ---- end of slicer settings"):
                    break
    if work_dir and os.path.exists(os.path.join(work_dir, "params_used.json")):
        with open(os.path.join(work_dir, "params_used.json"), encoding="utf-8") as fh:
            return np.asarray(json.load(fh)["deform"].get("PART_OFFSET", [0.0, 0.0, 0.0]), float)
    return np.zeros(3)


def _model(stl_path, part_offset=(0.0, 0.0, 0.0)):
    m = o3d.io.read_triangle_mesh(stl_path)
    V = np.asarray(m.vertices)
    shift = np.array([(V[:, 0].min() + V[:, 0].max()) / 2, (V[:, 1].min() + V[:, 1].max()) / 2, V[:, 2].min()])
    m.vertices = o3d.utility.Vector3dVector(V - (shift + np.asarray(part_offset, float)))
    return m


def surface_centroid(m):
    """Area-weighted centroid of the mesh's triangles (exact; what a uniform surface sample's mean estimates)."""
    V = np.asarray(m.vertices)
    T = V[np.asarray(m.triangles)]
    area = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
    return (area[:, None] * T.mean(1)).sum(0) / area.sum()


def measure(model, stl_path, mid, seg_len, rough_deg, part_offset=(0.0, 0.0, 0.0)):
    """mid, seg_len, rough_deg: per bead segment (quality.path_roughness's mid and rough, plus lengths).
    part_offset: the run's PART_OFFSET (part_offset_of)."""
    m = _model(stl_path, part_offset)
    scene = o3d.t.geometry.RaycastingScene()
    scene.add_triangles(o3d.t.geometry.TriangleMesh.from_legacy(m))
    q = o3d.core.Tensor(mid.astype(np.float32))
    # the sign by a 5-ray inside/outside vote: compute_signed_distance casts one ray, and on a mesh that isn't
    # watertight (Squirtle's STL is open and self-intersecting) it flips the sign of a few points deep inside
    sd = scene.compute_distance(q).numpy() * np.where(scene.compute_occupancy(q, nsamples=5).numpy() > 0.5, -1.0, 1.0)
    out = {"outside_mm": round(float(seg_len[sd > 1.0].sum()), 1)}
    near = np.abs(sd) < 0.7
    out["surface_rough_pct"] = round(100.0 * float((rough_deg[near] > 30).mean()), 2) if near.any() else None
    if model in SYMMETRY:
        a = np.radians(SYMMETRY[model])
        n = np.array([np.cos(a), np.sin(a), 0.0])
        c = surface_centroid(m)
        P = mid[::4]
        Pm = P - 2 * ((P - c) @ n)[:, None] * n
        out["asym_mm"] = round(float(np.percentile(cKDTree(mid).query(Pm)[0], 90)), 3)
    else:
        out["asym_mm"] = None
    return out
