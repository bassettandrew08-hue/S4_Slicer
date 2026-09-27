"""
Print-order support check: find extrusion that would be laid down in mid-air.

Works on the real-space toolpath (planar G-code mapped back through the deformation, exactly as the
mapper does it), in print order, one planar layer at a time. A point is "floating" if nothing printed in
an earlier layer (or the bed) is within `radius` mm of it. Every point also gets checked against the real
part, to confirm the toolpath stays inside the model.
"""
import re
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

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
          retraction_length=1.0):
    inp = meshio_s4.load_and_tetrahedralize(model_path, part_offset)
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

    surf = inp.extract_surface(algorithm="dataset_surface")
    inside = np.asarray(inp.find_containing_cell(real[ext])) >= 0
    _, cp = surf.find_closest_cell(real[ext], return_closest_point=True)
    outside_mm = np.where(inside, 0.0, np.linalg.norm(real[ext] - cp, axis=1))

    fl = np.nonzero(floating)[0]

    # Group floating points into runs along the toolpath and classify them by what the ends connect to:
    #   bridge     - supported extrusion on both ends (normal FDM bridging; sparse gyroid does this even in flat prints)
    #   cantilever - supported on one end only
    #   island     - not connected to supported extrusion at all: this plastic has nothing to stick to
    def supported(i):
        return 0 <= i < len(P) and ext[i] and not floating[i]
    runs = np.split(fl, np.nonzero(np.diff(fl) > 1)[0] + 1) if len(fl) else []
    classes = {"bridge": [0, 0.0, Counter()], "cantilever": [0, 0.0, Counter()], "island": [0, 0.0, Counter()]}
    for r in runs:
        a, b = supported(r[0] - 1), supported(r[-1] + 1)
        k = "bridge" if a and b else ("cantilever" if a or b else "island")
        chain = np.r_[r[0] - 1, r] if r[0] > 0 else r
        classes[k][0] += 1
        classes[k][1] += float(np.sum(np.linalg.norm(np.diff(real[chain], axis=0), axis=1)))
        classes[k][2].update(types[i] for i in r)

    return {
        "runs": {k: {"count": v[0], "length_mm": v[1], "types": v[2].most_common()} for k, v in classes.items()},
        "extruding_points": int(ext.sum()),
        "floating_points": int(len(fl)),
        "floating_pct": 100.0 * len(fl) / max(int(ext.sum()), 1),
        "floating_types": Counter(types[i] for i in fl).most_common(),
        "floating_gap_mm": (float(np.median(gap[fl])), float(gap[fl].max())) if len(fl) else (0.0, 0.0),
        "floating_layers": (int(layer[fl].min()), int(layer[fl].max())) if len(fl) else None,
        "floating_real_z": (float(real[fl, 2].min()), float(real[fl, 2].max())) if len(fl) else None,
        "outside_part_points": int((outside_mm > 0.5).sum()),
        "radius_mm": radius,
    }


def format_report(r):
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
