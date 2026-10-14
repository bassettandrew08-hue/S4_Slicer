"""
How smoothly the layer lines run over a model's visible outer surface (the "choppy hull" of the benchy).

Every printed layer is a level set of the deformed height: a layer line on the outer surface is a contour line of
h = deformed z on the undeformed surface, one every layer height. This measures those contour lines directly on the
tet mesh's surface triangles (h is linear on each, as in the mapping), in a region of the visible exterior:

  layer_turns       total turning of the layer lines in the region, in full turns (360 deg): the sum over surface
                    edges of (layer lines crossing the edge) x (how much their direction changes there, measured in
                    the unfolded surface). Smooth layers on a smooth surface turn little; a closed ring counts 1;
                    lines that zig-zag at tet faces or wander where layers graze the surface count a lot.
  layer_kink_turns  the same, only for direction changes over 15 deg at a single edge (kinks, not curves)

The planar print of the undeformed model (h = z) is printed as `flat` for scale. Regions are visible-exterior
surface (orthographic rays from +-x, +-y and +z) inside a box, in the part's frame (centred in x and y, z = 0 at its
bottom, as meshio_s4.load_and_tetrahedralize places it).

    venv\\Scripts\\python tools\\layer_lines.py build/suite/base [build/suite/try ...] [--model NAME]

Each argument is a model_suite --out directory (or a run's work directory holding deformed_points.npy and
params_used.json).
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402

# region name: x, y, z box of the visible exterior surface
REGIONS = {"benchy upsidedown tilted": {"hull": ((-40, 40), (-40, 40), (26.5, 60))}}

_MESH = {}


def _surface(stl, part_offset, simplify_error, max_volume, box):
    """The tet mesh's surface (vertices, faces, original point ids, unit normals, areas) and the faces of the
    visible exterior inside box, plus the interior edges between two such faces on smooth surface."""
    key = (stl, tuple(part_offset), simplify_error, max_volume, box)
    if key in _MESH:
        return _MESH[key]
    import open3d as o3d
    from s4 import meshio_s4
    grid = meshio_s4.load_and_tetrahedralize(stl, part_offset, simplify_error, max_volume)
    surf = grid.extract_surface(algorithm="dataset_surface")
    F = surf.faces.reshape(-1, 4)[:, 1:].astype(np.int64)
    orig = np.asarray(surf.point_data["vtkOriginalPointIds"])
    V = np.asarray(surf.points)
    nf = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    area = np.linalg.norm(nf, axis=1) / 2
    nf /= np.maximum(2 * area[:, None], 1e-300)
    scene = o3d.t.geometry.RaycastingScene()
    scene.add_triangles(o3d.core.Tensor(V.astype(np.float32)), o3d.core.Tensor(F.astype(np.uint32)))
    lo, hi = V.min(0) - 1, V.max(0) + 1
    vis = np.zeros(len(F), bool)
    step = 0.1
    for ax in range(3):  # rays along -ax and +ax (z: from above only)
        u, w = [k for k in range(3) if k != ax]
        A, B = np.meshgrid(np.arange(lo[u], hi[u], step), np.arange(lo[w], hi[w], step))
        for sign in ((1, -1) if ax < 2 else (1,)):
            rays = np.zeros((A.size, 6), np.float32)
            rays[:, u], rays[:, w] = A.ravel(), B.ravel()
            rays[:, ax] = hi[ax] + 10 if sign > 0 else lo[ax] - 10
            rays[:, 3 + ax] = -sign
            ans = scene.cast_rays(o3d.core.Tensor(rays))
            hit = ans["t_hit"].numpy() < 1e30
            vis[ans["primitive_ids"].numpy().astype(np.int64)[hit]] = True
    ctr = V[F].mean(1)
    region = vis.copy()
    for k, (a, b) in enumerate(box):
        region &= (ctr[:, k] >= a) & (ctr[:, k] <= b)
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    fo = np.tile(np.arange(len(F)), 3)
    e = np.sort(e, axis=1)
    o = np.lexsort((e[:, 1], e[:, 0]))
    e, fo = e[o], fo[o]
    i = np.nonzero(np.all(e[1:] == e[:-1], axis=1))[0]
    E, f1, f2 = e[i], fo[i], fo[i + 1]
    # both faces in the region, and the surface smooth there (no model edges: their layer lines turn anyway)
    m = region[f1] & region[f2] & (np.einsum("ij,ij->i", nf[f1], nf[f2]) > np.cos(np.radians(30)))
    out = dict(F=F, orig=orig, V=V, nf=nf, area=area, region=region, E=E[m], f1=f1[m], f2=f2[m])
    _MESH[key] = out
    return out


def turns(s, h, layer_height):
    """(layer_turns, layer_kink_turns) of the contour lines of h (per surface vertex) every layer_height."""
    F, V, nf, area, E, f1, f2 = s["F"], s["V"], s["nf"], s["area"], s["E"], s["f1"], s["f2"]
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    h0, h1, h2 = h[F[:, 0]], h[F[:, 1]], h[F[:, 2]]
    grad = ((h1 - h0)[:, None] * np.cross(nf, p0 - p2) + (h2 - h0)[:, None] * np.cross(nf, p1 - p0)) \
        / np.maximum(2 * area[:, None], 1e-300)
    ed = V[E[:, 1]] - V[E[:, 0]]
    ed /= np.linalg.norm(ed, axis=1)[:, None]

    def frame(f):  # contour direction n x grad in the (edge, n x edge) frame: the same frame on both faces after
        t = np.cross(nf[f], grad[f])  # unfolding the hinge
        return np.einsum("ij,ij->i", t, ed) + 1j * np.einsum("ij,ij->i", t, np.cross(nf[f], ed))
    ang = np.abs(np.angle(frame(f1) * np.conj(frame(f2))))
    t = np.abs(h[E[:, 0]] - h[E[:, 1]]) / layer_height * ang
    return float(t.sum() / (2 * np.pi)), float(t[ang > np.radians(15)].sum() / (2 * np.pi))


def _layer_height(planar):
    if planar and os.path.exists(planar):
        with open(planar, errors="replace") as fh:
            for line in fh:
                m = re.match(r";Layer height: ([\d.]+)", line)
                if m:
                    return float(m.group(1))
                if line.startswith(";LAYER:"):
                    break
    return None


def measure(model, work, layer_height=None, flat=False):
    """{f"layer_turns_{region}": .., f"layer_kink_turns_{region}": ..} for a run's work directory (holding
    deformed_points.npy and params_used.json). flat=True: the planar print of the undeformed model instead."""
    if model not in REGIONS:
        return {}
    from s4.params import expand_iterations
    prof = json.load(open(os.path.join(work, "params_used.json")))
    off, iters = expand_iterations(prof["deform"])
    lh = layer_height or _layer_height(os.path.join(work, f"{model}_deformed_tet.gcode")) or 0.13
    Pd = None if flat else np.load(os.path.join(work, "deformed_points.npy"))
    res = {}
    for name, box in REGIONS[model].items():
        s = _surface(os.path.join(HERE, "input_models", f"{model}.stl"), tuple(off), iters[0]["SURFACE_SIMPLIFY_ERROR"],
                     iters[0].get("TET_MAX_VOLUME", 0.0), box)
        h = s["V"][:, 2] if flat else Pd[s["orig"], 2]
        a, b = turns(s, h, lh)
        res[f"layer_turns_{name}"], res[f"layer_kink_turns_{name}"] = round(a, 1), round(b, 1)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1], formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("runs", nargs="+", help="model_suite --out directories or run work directories")
    ap.add_argument("--model", default="benchy upsidedown tilted")
    a = ap.parse_args()
    first = True
    for r in a.runs:
        work = os.path.join(r, a.model) if os.path.isdir(os.path.join(r, a.model)) else r
        if first:
            print(f"{'flat (planar print of the undeformed model)':60s} {measure(a.model, work, flat=True)}")
            first = False
        print(f"{r:60s} {measure(a.model, work)}")


if __name__ == "__main__":
    main()
