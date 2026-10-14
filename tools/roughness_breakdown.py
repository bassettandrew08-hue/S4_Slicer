"""
Where does the path roughness (s4/quality.py path_roughness, the suite's "zigzag") come from?

Re-maps a finished run's planar G-code with the run's deformed mesh, keeps track of which planar point every output
move came from, and sorts every rough bead segment (> 30 deg, neighbour also sharp) by cause:

  cura_rough     the planar (flat) path is itself rough there: Cura's own pattern (skin/infill zig-zags, line ends,
                 U-turns). It turns just as sharply in a flat print
  cura_corner    the planar path has one sharp corner there (a real corner) and the mapping adds the second turn
  in_move        the turn is inside one output move: the sim interpolates C/X/Z/B linearly, so a move whose C or B
                 changes a lot is drawn as a curve (near the turntable axis, or with fast B changes)
  rot_split      a rotation-split step (ROTATION_MAX_DELTA)
  travel         next to a re-entry / plunge / retract line
  tet_face       the planar path is straight, the mapped one kinks: the move crosses a tet face (the deformation is
                 linear per tet, so straight lines bend at faces)
  same_tet       straight planar path inside one tet that still turns (should be rare)

    venv\\Scripts\\python tools\\roughness_breakdown.py build/suite/base "pi 3mm" [--set KEY=VALUE ...]

The run's directory must hold <model>/<model>_deformed_tet.gcode and deformed_points.npy (a model_suite --out).
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402

from s4 import profile as profiles, fast_map, meshio_s4  # noqa: E402
from s4.params import expand_iterations  # noqa: E402
from s4.pipeline import planar_retraction  # noqa: E402
from s4.geometry import nozzle_tip  # noqa: E402

TYPES = ["WALL-OUTER", "WALL-INNER", "SKIN", "FILL", "OTHER"]


def planar_types(path, seg):
    """Feature type (index into TYPES) per segmented planar point, same segmentation as fast_map.read_gcode_points."""
    word = fast_map._WORD
    pos = np.array([0., 0., 20.])
    out = []
    t = 4
    for raw in open(path).read().splitlines():
        if raw.startswith(";TYPE:"):
            name = raw[6:].strip()
            t = TYPES.index(name) if name in TYPES else 4
            continue
        code = raw.split(";", 1)[0]
        if "(" in code:
            code = re.sub(r"\([^)]*\)", "", code)
        words = word.findall(code)
        motion = None
        new = pos.copy()
        for L, num in words:
            L = L.upper()
            if L == "G" and float(num) in (0, 1):
                motion = 1
            elif L in "XYZ":
                new["XYZ".index(L)] = float(num)
        if motion is None:
            continue
        d = np.linalg.norm(new - pos)
        n = int(-(-d // seg)) if d > 0 else 1
        out.extend([t] * n)
        pos = new
    return np.array(out)


def turn(d0, d1):
    L0, L1 = np.linalg.norm(d0, axis=-1), np.linalg.norm(d1, axis=-1)
    c = np.einsum("...i,...i", d0, d1) / np.maximum(L0 * L1, 1e-12)
    return np.where((L0 > 1e-6) & (L1 > 1e-6), np.degrees(np.arccos(np.clip(c, -1, 1))), 0.0)


def planar_rough_points(path, seg=0.6, ret=1.0):
    """Mid points of the rough segments of a planar G-code (a flat print scored the same way)."""
    g = fast_map.read_gcode_points(path, seg)
    P = g["position"]
    ext = np.array([e if e is not None else 0.0 for e in g["extrusion"]])
    prn = np.array([c == "G01" for c in g["command"]]) & (ext > 0) & ~np.isclose(ext, ret)
    d = np.diff(P, axis=0, prepend=P[:1])
    raw = np.zeros(len(P))
    raw[1:] = np.where(prn[1:] & prn[:-1], turn(d[:-1], d[1:]), 0.0)
    nb = np.maximum(np.r_[0, raw[:-1]], np.r_[raw[1:], 0])
    mid = P - d / 2
    return mid[np.minimum(raw, nb) > 30]


def waviness(gcode_path, planar_path, src, seg=0.6, ret=1.0, nozzle_offset=42.0, threshold=10.0):
    """Joints of the final path that turn more than `threshold` deg where the planar path runs straight on (planar
    turn < 2 deg): what the mapping adds as wavy walls. src: the planar point of every G00/G01 line the mapper wrote
    (the run's map_src.npy). Returns (count, straight joints)."""
    W = re.compile(r"([CXZBE])(-?\d+(?:\.\d*)?|-?\.\d+)")
    C = X = Z = B = 0.0
    rows, owner = [], []
    m = -1
    with open(gcode_path, errors="replace") as fh:
        for line in fh:
            if not line.startswith(("G0", "G1")):
                continue
            mapper = line.startswith(("G00 ", "G01 "))
            m += mapper
            w = dict(W.findall(line.split(";", 1)[0]))
            c, x, z, b = (float(w.get(k, v)) for k, v in (("C", C), ("X", X), ("Z", Z), ("B", B)))
            if mapper and line.startswith("G01") and float(w.get("E", 0)) > 1e-6 and (c, x, z, b) != (C, X, Z, B):
                rows.append((C, X, Z, B, c, x, z, b))
                owner.append(m)
            C, X, Z, B = c, x, z, b
    if m + 1 != len(src):
        raise ValueError(f"map_src has {len(src)} moves, the G-code {m + 1}")
    a = np.array(rows)
    A0, A1 = a[:, :4], a[:, 4:]
    k = np.clip(np.maximum(np.ceil(np.abs(A1[:, 0] - A0[:, 0]) / 3), np.ceil(np.abs(A1[:, 3] - A0[:, 3]) / 6)), 1, 32)
    k = k.astype(int)
    sm = np.repeat(np.arange(len(a)), k)
    jj = np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)

    def tip(f):
        s = A0[sm] + (A1[sm] - A0[sm]) * f[:, None]
        return np.c_[nozzle_tip(s[:, 0], s[:, 1], s[:, 2], s[:, 3], nozzle_offset)]
    p0, p1 = tip(jj / k[sm]), tip((jj + 1) / k[sm])
    dd = p1 - p0
    Ls = np.linalg.norm(dd, axis=1)
    raw = np.zeros(len(p0))
    cont = (np.linalg.norm(p0[1:] - p1[:-1], axis=1) <= 0.05) & (Ls[1:] > 1e-6) & (Ls[:-1] > 1e-6)
    raw[1:] = np.where(cont, turn(dd[:-1], dd[1:]), 0.0)
    g = fast_map.read_gcode_points(planar_path, seg)
    P = g["position"]
    ext = np.array([e if e is not None else 0.0 for e in g["extrusion"]])
    prn = np.array([c == "G01" for c in g["command"]]) & (ext > 0) & ~np.isclose(ext, ret)
    d = np.diff(P, axis=0, prepend=P[:1])
    praw = np.zeros(len(P))
    praw[1:] = np.where(prn[1:] & prn[:-1], turn(d[:-1], d[1:]), 0.0)
    fq = np.nonzero(jj == 0)[0]
    si = np.asarray(src)[np.asarray(owner)[sm[fq]]]
    straight = (praw[si] < 2.0) & prn[si] & (raw[fq] > 0)
    return int((raw[fq][straight] > threshold).sum()), int(straight.sum())


def analyse(run_dir, model, sets=(), mp_over=None, verbose=True, flat=None):
    stl = os.path.join(HERE, "input_models", f"{model}.stl")
    work = os.path.join(run_dir, model)
    planar = os.path.join(work, f"{model}_deformed_tet.gcode")
    prof, _ = profiles.resolve(stl, sets=sets)
    p = profiles.deform_params(prof)
    off, iters = expand_iterations(p)
    grid = meshio_s4.load_and_tetrahedralize(stl, off, iters[0]["SURFACE_SIMPLIFY_ERROR"], iters[0].get("TET_MAX_VOLUME"))
    cells = meshio_s4.cells_of(grid)
    in_pts = np.asarray(grid.points)
    in_c = np.asarray(grid.cell_centers().points)
    dp = np.load(os.path.join(work, "deformed_points.npy"))
    dgrid = meshio_s4.deformed_grid_from_points(grid, dp)
    dc = np.asarray(dgrid.cell_centers().points)
    mp = dict(prof["map"])
    mp.update(mp_over or {})
    if mp["RETRACTION_LENGTH"] is None:
        mp["RETRACTION_LENGTH"] = planar_retraction(planar)
    store = {}
    finder = meshio_s4.make_find_cells(dgrid)

    def find(P):
        store["c"] = finder(P)
        return store["c"]
    pts, _ = fast_map.map_gcode(cells, in_pts, in_c, dp, dc, find, planar, mp)
    containing = store["c"]
    g = fast_map.read_gcode_points(planar, mp["SEG_SIZE"])
    P = g["position"]
    ptype = planar_types(planar, mp["SEG_SIZE"])
    assert len(ptype) == len(P), (len(ptype), len(P))
    RET = mp["RETRACTION_LENGTH"]
    # planar turn per planar segment (segment i = P[i-1] -> P[i]), printing segments only
    ext = np.array([e if e is not None else 0.0 for e in g["extrusion"]])
    prn = np.array([c == "G01" for c in g["command"]]) & (ext > 0) & ~np.isclose(ext, RET)
    d = np.diff(P, axis=0, prepend=P[:1])
    praw = np.zeros(len(P))
    praw[1:] = np.where(prn[1:] & prn[:-1], turn(d[:-1], d[1:]), 0.0)
    pnb = np.maximum(np.r_[0, praw[:-1]], np.r_[praw[1:], 0])
    prough = np.minimum(praw, pnb)

    # output moves -> machine coordinates exactly as write_gcode, then the sim's bead segments
    L = mp["NOZZLE_OFFSET"]
    O = np.asarray(pts["position"], float)
    rot = np.asarray(pts["rotation"], float)
    keep = np.nonzero(~np.isnan(O).any(1) & ~(O[:, 2] < 0))[0]
    O, rot = O[keep], rot[keep]
    zh = np.asarray(pts["travelling"], float)[keep]
    r = np.linalg.norm(O[:, :2], axis=1) - np.sin(rot) * (L + zh)
    th = np.arctan2(O[:, 1], O[:, 0])
    z = O[:, 2] + (np.cos(rot) - 1) * (L + zh) + zh
    dl = np.diff(np.r_[0.0, th])
    dl = np.where(dl > np.pi, dl - 2 * np.pi, dl)
    dl = np.where(dl < -np.pi, dl + 2 * np.pi, dl)
    M = np.c_[np.round(np.rad2deg(np.cumsum(dl)), 5), np.round(r, 5), np.round(z, 5), np.round(np.rad2deg(rot), 5)]
    cmd = [pts["command"][j] for j in keep]
    e = [pts["extrusion"][j] for j in keep]
    src = np.asarray(pts["src"])[keep]
    kind = np.asarray(pts["kind"])[keep]
    prev = np.r_[[[0, 0, 20, 0]], M[:-1]]
    extr = np.array([c == "G01" and x is not None and round(x, 4) > 1e-6 for c, x in zip(cmd, e)])
    extr &= np.any(M != prev, axis=1)
    mv = np.nonzero(extr)[0]
    A0, A1 = prev[mv], M[mv]
    k = np.clip(np.maximum(np.ceil(np.abs(A1[:, 0] - A0[:, 0]) / 3), np.ceil(np.abs(A1[:, 3] - A0[:, 3]) / 6)), 1, 32)
    k = k.astype(int)
    sm = np.repeat(np.arange(len(mv)), k)
    jj = np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)

    def tip(f):
        s = A0[sm] + (A1[sm] - A0[sm]) * f[:, None]
        return np.c_[nozzle_tip(s[:, 0], s[:, 1], s[:, 2], s[:, 3], L)]
    p0, p1 = tip(jj / k[sm]), tip((jj + 1) / k[sm])
    dd = p1 - p0
    Ls = np.linalg.norm(dd, axis=1)
    raw = np.zeros(len(p0))
    cont = (np.linalg.norm(p0[1:] - p1[:-1], axis=1) <= 0.05) & (Ls[1:] > 1e-6) & (Ls[:-1] > 1e-6)
    raw[1:] = np.where(cont, turn(dd[:-1], dd[1:]), 0.0)
    nb = np.maximum(np.r_[0, raw[:-1]], np.r_[raw[1:], 0])
    rough = np.minimum(raw, nb)
    R = np.nonzero(rough > 30)[0]

    # per rough sub-segment: which move (output index), which planar point
    om = mv[sm[R]]          # index into keep-arrays
    first = jj[R] == 0
    i = src[om]
    n = len(P)
    win = lambda a: np.maximum.reduce([a[np.clip(i + o, 0, n - 1)] for o in (-1, 0, 1)])
    cls = np.full(len(R), "", dtype=object)
    near_kind = np.maximum.reduce([np.r_[kind[1:], 0], kind, np.r_[0, kind[:-1]]])
    cls[(win(prough) > 30)] = "cura_rough"
    m = cls == ""
    cls[m & (win(praw) > 30)] = "cura_corner"
    m = cls == ""
    # the turn of this sub-segment sits inside its move, or the neighbour that makes it rough does
    cls[m & ~first] = "in_move"
    m = cls == ""
    cls[m & (kind[om] == 1)] = "rot_split"
    m = cls == ""
    cls[m & np.isin(near_kind[om], (2, 3))] = "travel"
    m = cls == ""
    cprev = containing[np.clip(i - 1, 0, n - 1)]
    cnext = containing[np.clip(i + 1, 0, n - 1)]
    face = (cprev != containing[i]) | (cnext != containing[i])
    cls[m & face] = "tet_face"
    cls[cls == ""] = "same_tet"
    # surface rough share, as tools/shape_metrics
    from tools import shape_metrics
    import open3d as o3d
    mesh = shape_metrics._model(stl)
    scene = o3d.t.geometry.RaycastingScene()
    scene.add_triangles(o3d.t.geometry.TriangleMesh.from_legacy(mesh))
    mid = (p0 + p1) / 2
    sd = scene.compute_signed_distance(o3d.core.Tensor(mid.astype(np.float32))).numpy()
    near = np.abs(sd) < 0.7
    res = {"segments": len(p0), "count": len(R), "surface_rough_pct": round(100 * float((rough[near] > 30).mean()), 2),
           "by_cause": {}, "by_type": {}, "surface_by_cause": {}}
    t = ptype[i]
    for c in ("cura_rough", "cura_corner", "in_move", "rot_split", "travel", "tet_face", "same_tet"):
        res["by_cause"][c] = int((cls == c).sum())
        res["surface_by_cause"][c] = int(((cls == c) & near[R]).sum())
    for ti, tn in enumerate(TYPES):
        res["by_type"][tn] = {c: int(((cls == c) & (t == ti)).sum()) for c in res["by_cause"] if
                              ((cls == c) & (t == ti)).any()}
    res["surface_near_segments"] = int(near.sum())
    res["smooth"] = dict(fast_map.SMOOTH_STATS) if mp.get("MAP_SMOOTH") else None
    # mapping wobble: turns of the final path where the planar path goes straight on (planar turn < 2 deg), at the
    # joints between output moves. Below the 30 deg roughness threshold, but what shows as a wavy wall
    fq = np.nonzero(jj == 0)[0]
    si = src[mv[sm[fq]]]
    straight = (praw[si] < 2.0) & prn[si] & (raw[fq] > 0)
    tw = raw[fq][straight]
    res["wobble"] = {"joints": int(straight.sum()), "mean_deg": round(float(tw.mean()), 3) if len(tw) else 0,
                     "gt5": int((tw > 5).sum()), "gt10": int((tw > 10).sum()), "gt20": int((tw > 20).sum())}
    rr = np.linalg.norm(p0[fq][straight][:, :2], axis=1)
    res["wobble"]["gt10_r_lt5mm"] = int(((tw > 10) & (rr < 5)).sum())
    dC = np.abs(A1[:, 0] - A0[:, 0])[sm[fq][straight]]
    res["wobble"]["gt10_dC_gt3deg"] = int(((tw > 10) & (dC > 3)).sum())
    if flat and os.path.exists(flat):
        # is the same spot rough in a flat print of the undeformed model (the model's own geometry)?
        from scipy.spatial import cKDTree
        fr = planar_rough_points(flat, mp["SEG_SIZE"], RET)
        dist = cKDTree(fr).query(mid[R])[0] if len(fr) else np.full(len(R), np.inf)
        inflat = dist < 1.0
        res["flat_rough_count"] = len(fr)
        res["also_rough_in_flat"] = {c: int(((cls == c) & inflat).sum()) for c in res["by_cause"]}
    res["planar_rough_count"] = int((prough > 30).sum())
    if verbose:
        print(model, {k: v for k, v in res.items() if k not in ("by_type",)})
        for tn, v in res["by_type"].items():
            print("   ", tn, v)
    res["_arrays"] = dict(R=R, cls=cls, mid=mid, rough=rough, near=near, sd=sd, i=i, om=om)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("models", nargs="+")
    ap.add_argument("--set", action="append", default=[])
    a = ap.parse_args()
    for mdl in a.models:
        analyse(a.run_dir, mdl, a.set, flat=os.path.join(HERE, "build", "flat", f"{mdl}_flat.gcode"))


if __name__ == "__main__":
    main()
