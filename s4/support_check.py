"""
Print-order support check: find extrusion that would be laid down in mid-air.
(The poles / along-the-nozzle-axis check of the final 4-axis G-code is in s4/quality.py.)

Works on the real-space toolpath (planar G-code mapped back through the deformation, exactly as the
mapper does it), in print order, one planar layer at a time. A point is "floating" if nothing printed in
an earlier layer (or the bed) is within `radius` mm of it. Every point also gets checked against the real
part, to confirm the toolpath stays inside the model.
"""
import re
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from . import fast_map, meshio_s4

_WORD = re.compile(r"([A-Za-z])\s*(-?(?:\d+\.?\d*|\.\d+))")


def _planar_layers(planar_path, seg_size):
    """Planar layer number and ;TYPE: for every segmented point read_gcode_points produces."""
    layer_ids, types = [], []
    layer, typ = -1, "?"
    pos = np.array([0., 0., 20.])
    with open(planar_path) as fh:
        for raw in fh:
            if raw.startswith(";LAYER:"):
                layer = int(raw[7:])
            elif raw.startswith(";TYPE:"):
                typ = raw[6:].strip()
            words = _WORD.findall(raw.split(";", 1)[0])
            if not words:
                continue
            g = [float(n) for L, n in words if L.upper() == "G"]
            if not g or g[0] not in (0.0, 1.0):
                continue
            prev = pos.copy()
            for L, n in words:
                k = "XYZ".find(L.upper())
                if k >= 0:
                    pos[k] = float(n)
            dist = np.linalg.norm(pos - prev)
            n = int(-(-dist // seg_size)) if dist > 0 else 1
            layer_ids += [layer] * n
            types += [typ] * n
    return np.array(layer_ids), types


def check(model_path, deformed_points, planar_path, part_offset=(0., 0., 0.), radius=1.0, seg_size=0.6,
          retraction_length=1.0, input_grid=None):
    """input_grid: undeformed tet mesh to use instead of re-running tetgen on model_path."""
    inp = input_grid if input_grid is not None else meshio_s4.load_and_tetrahedralize(model_path, part_offset)
    cells = meshio_s4.cells_of(inp)
    P0 = np.asarray(inp.points)
    Pd = np.asarray(deformed_points)
    d = meshio_s4.deformed_grid_from_points(inp, Pd)

    g = fast_map.read_gcode_points(planar_path, seg_size)
    P = g["position"]
    layer, types = _planar_layers(planar_path, seg_size)
    assert len(layer) == len(P)
    cont = meshio_s4.make_find_cells(d)(P)
    cv = Pd[cells[cont]]
    total = fast_map.tet_volumes(*[cv[:, i] for i in range(4)])
    bary = np.stack([fast_map.tet_volumes(P, cv[:, 1], cv[:, 2], cv[:, 3]),
                     fast_map.tet_volumes(P, cv[:, 0], cv[:, 2], cv[:, 3]),
                     fast_map.tet_volumes(P, cv[:, 0], cv[:, 1], cv[:, 3]),
                     fast_map.tet_volumes(P, cv[:, 0], cv[:, 1], cv[:, 2])], axis=1) / total[:, None]
    real = P - np.sum((Pd - P0)[cells[cont]] * bary[:, :, None], axis=1)
    ext = np.array([e is not None and e > 0 and abs(e) != retraction_length for e in g["extrusion"]])
    ext &= ~(bary.sum(axis=1) > 1.01)  # the mapper drops these

    surf = inp.extract_surface(algorithm="dataset_surface")
    inside = np.asarray(inp.find_containing_cell(real[ext])) >= 0
    _, cp = surf.find_closest_cell(real[ext], return_closest_point=True)
    outside_mm = np.where(inside, 0.0, np.linalg.norm(real[ext] - cp, axis=1))
    r = analyse(real, ext, layer, types, radius)
    r["outside_part_points"] = int((outside_mm > 0.5).sum())
    return r


def check_planar(planar_path, radius=1.0, seg_size=0.6, retraction_length=1.0):
    """Debugging helper, not called by the pipeline: the same check on Cura's planar toolpath in the deformed shape
    (no mapping back), to see whether an island is born in the deformation or in the mapping."""
    g = fast_map.read_gcode_points(planar_path, seg_size)
    layer, types = _planar_layers(planar_path, seg_size)
    ext = np.array([e is not None and e > 0 and abs(e) != retraction_length for e in g["extrusion"]])
    r = analyse(g["position"], ext, layer, types, radius)
    r["outside_part_points"] = 0
    return r


def analyse(real, ext, layer, types, radius):
    """Print-order support analysis of an extrusion path (points `real`, in print order)."""
    P = real
    # supported = the bed is within radius, or some point from an earlier planar layer is
    ei = np.nonzero(ext)[0]
    Q = real[ei]
    Ly = layer[ei]
    tree = cKDTree(Q)
    unsupported = Q[:, 2] > radius
    cand = np.nonzero(unsupported)[0]
    for k, nbrs in zip(cand, tree.query_ball_point(Q[cand], radius)):
        if len(nbrs) and Ly[np.asarray(nbrs)].min() < Ly[k]:
            unsupported[k] = False
    floating = np.zeros(len(P), bool)
    floating[ei[unsupported]] = True
    gap = np.full(len(P), np.nan)
    for k in np.nonzero(unsupported)[0]:  # exact gap only for the (few) floating points
        earlier = Q[Ly < Ly[k]]
        g = np.min(np.linalg.norm(earlier - Q[k], axis=1)) if len(earlier) else np.inf
        gap[ei[k]] = min(g, Q[k, 2])

    fl = np.nonzero(floating)[0]
    ug_mm, regions, ug_len = grounded_regions(real, ext, layer, types, radius)

    # Group floating points into runs along the toolpath and classify them by what the ends connect to:
    #   bridge     - supported extrusion on both ends (normal FDM bridging; sparse gyroid does this even in flat prints)
    #   cantilever - supported on one end only
    #   island     - not connected to supported extrusion at all: this plastic has nothing to stick to
    def supported(i):
        return 0 <= i < len(P) and ext[i] and not floating[i]
    runs = np.split(fl, np.nonzero(np.diff(fl) > 1)[0] + 1) if len(fl) else []
    classes = {"bridge": [0, 0.0, Counter(), 0.0], "cantilever": [0, 0.0, Counter(), 0.0], "island": [0, 0.0, Counter(), 0.0]}
    for r in runs:
        a, b = supported(r[0] - 1), supported(r[-1] + 1)
        k = "bridge" if a and b else ("cantilever" if a or b else "island")
        chain = np.r_[r[0] - 1, r] if r[0] > 0 else r
        length = float(np.sum(np.linalg.norm(np.diff(real[chain], axis=0), axis=1)))
        classes[k][0] += 1
        classes[k][1] += length
        classes[k][2].update(types[i] for i in r)
        # walls and skin share of the run: sparse infill floats a little even in a flat print (a flat cube with
        # this profile shows FILL cantilevers too), walls and skin don't
        classes[k][3] += length * float(np.mean([types[i] in STRUCTURAL for i in r]))

    return {
        "ungrounded_mm": ug_mm,
        "ungrounded_real_mm": ug_len["real_mm"],  # true path length
        "ungrounded_wall_mm": ug_len["wall_mm"],  # of it, walls and skin (not sparse infill)
        "regions": regions,
        "runs": {k: {"count": v[0], "length_mm": v[1], "wall_mm": v[3], "types": v[2].most_common()}
                 for k, v in classes.items()},
        "extruding_points": int(ext.sum()),
        "floating_points": int(len(fl)),
        "floating_pct": 100.0 * len(fl) / max(int(ext.sum()), 1),
        "floating_types": Counter(types[i] for i in fl).most_common(),
        "floating_gap_mm": (float(np.median(gap[fl])), float(gap[fl].max())) if len(fl) else (0.0, 0.0),
        "floating_layers": (int(layer[fl].min()), int(layer[fl].max())) if len(fl) else None,
        "floating_real_z": (float(real[fl, 2].min()), float(real[fl, 2].max())) if len(fl) else None,
        "radius_mm": radius,
    }


STRUCTURAL = ("WALL-OUTER", "WALL-INNER", "SKIN")  # Cura ;TYPE: lines that carry the part (not FILL)


def grounded_regions(P, ext, layer, types, R=1.0, reach=5.0, seg_mm=0.3):
    """Extrusion with no chain of support down to the bed.

    A point is grounded if the bed is within R, if a grounded point from an EARLIER layer is within R, or if it
    lies on the same extrusion run within `reach` mm of such a point (bridging / overhang reach). Everything else is
    ungrounded: material printed in mid-air, plus anything later stacked on top of it. Ungrounded points are grouped
    into regions (connected within R). Returns (total ungrounded mm, regions sorted by size)."""
    ei = np.nonzero(ext & (layer >= 0))[0]
    if len(ei) == 0:
        return 0.0, [], {"real_mm": 0.0, "wall_mm": 0.0}
    Q = P[ei]; Ly = layer[ei]; n = len(ei)
    pairs = cKDTree(Q).query_pairs(R, output_type="ndarray")
    a, b = pairs[:, 0], pairs[:, 1]
    fwd = Ly[a] > Ly[b]; bwd = Ly[b] > Ly[a]
    A = csr_matrix((np.ones(int(fwd.sum() + bwd.sum()), np.int8), (np.r_[a[fwd], b[bwd]], np.r_[b[fwd], a[bwd]])), shape=(n, n))
    run_id = np.cumsum(np.r_[1, (np.diff(ei) != 1) | (np.diff(Ly) != 0)])
    seg = np.r_[0, np.linalg.norm(np.diff(Q, axis=0), axis=1)] * (np.r_[0, np.diff(run_id)] == 0)
    s_along = np.cumsum(seg)
    g = np.zeros(n, bool)
    order = np.argsort(Ly, kind="stable")
    bounds = np.r_[0, np.cumsum(np.bincount(Ly - Ly.min()))]
    for k in range(len(bounds) - 1):
        idx = order[bounds[k]:bounds[k + 1]]
        if len(idx) == 0:
            continue
        direct = (Q[idx, 2] < R) | (np.asarray(A[idx] @ g.astype(np.int8)).ravel() > 0)
        gl = direct.copy()
        rid = run_id[idx]; sa = s_along[idx]
        for r in np.unique(rid[direct]):
            m = rid == r
            dist = np.min(np.abs(sa[m][:, None] - sa[m & direct][None, :]), axis=1)
            gl[np.nonzero(m)[0]] |= dist <= reach
        g[idx] = gl
    ug = ~g
    # path length each point stands for: half of the segments on either side of it within its run (points are ~0.5 mm
    # apart, so the old "count x seg_mm" understated the length)
    # at least seg_mm per point, so plastic piled into one spot (a blob: points with no path length) still counts
    point_len = 0.5 * (seg + np.r_[seg[1:], 0.0])
    point_len[point_len < 0.02] = seg_mm  # plastic piled into one spot (a blob) still counts; real points are further apart
    typ = np.array([types[i] for i in ei], dtype=object)
    structural = np.isin(typ, STRUCTURAL)  # walls, skin: sparse infill floats a little even in flat prints
    real_mm = float(point_len[ug].sum()); wall_mm = float(point_len[ug & structural].sum())
    keep = ug[a] & ug[b]
    C = csr_matrix((np.ones(int(keep.sum())), (a[keep], b[keep])), shape=(n, n))
    _, lab = connected_components(C, directed=False)
    regions = []
    for c in np.unique(lab[ug]):
        m = ug & (lab == c)
        L = Ly[m]
        regions.append({"mm": float(point_len[m].sum()), "layers": (int(L.min()), int(L.max())),
                        "centre": [round(float(v), 1) for v in Q[m].mean(0)],
                        "types": Counter(types[i] for i in ei[m]).most_common(2)})
    regions.sort(key=lambda r: -r["mm"])
    return float(ug.sum() * seg_mm), regions, {"real_mm": real_mm, "wall_mm": wall_mm}


def format_report(r):
    regs = r.get("regions", [])
    s = (f"[support] ungrounded (no support chain to the bed): {r.get('ungrounded_real_mm', 0):.0f} mm of extrusion in "
         f"{len(regs)} regions, {r.get('ungrounded_wall_mm', 0):.0f} mm of it walls/skin (the rest sparse infill)")
    for g in regs[:3]:
        s += (f"\n[support]   ~{g['mm']:.0f} mm, layers {g['layers'][0]}-{g['layers'][1]}, at {g['centre']}, "
              f"{', '.join(f'{t} {n}' for t, n in g['types'])}")
    s += "\n" + _format_local(r)
    return s


def _format_local(r):
    s = (f"[support] {r['floating_points']} of {r['extruding_points']} extruded points ({r['floating_pct']:.2f}%) "
         f"are more than {r['radius_mm']} mm from anything printed earlier")
    if r["floating_points"]:
        s += (f" (real z {r['floating_real_z'][0]:.1f}-{r['floating_real_z'][1]:.1f} mm, "
              f"gap up to {r['floating_gap_mm'][1]:.1f} mm):")
        notes = {"bridge": "anchored both ends, normal bridging",
                 "cantilever": "anchored one end",
                 "island": "UNANCHORED, nothing to stick to"}
        for k in ("island", "cantilever", "bridge"):
            v = r["runs"][k]
            if v["count"]:
                s += (f"\n[support]   {k:10s} {v['count']:4d} runs, {v['length_mm']:6.1f} mm "
                      f"({notes[k]}; {', '.join(f'{t} {n}' for t, n in v['types'][:3])})")
    if r["outside_part_points"]:
        s += f"\n[support]   {r['outside_part_points']} points map >0.5 mm outside the model"
    return s

