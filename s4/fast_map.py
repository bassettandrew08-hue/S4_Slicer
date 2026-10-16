"""
Fast G-code mapping (notebook cells 16-18), numerically equivalent to s4/reference.py.

What changed vs. the notebook (all output-preserving):
  * cell 16's calculate_tet_attributes(deformed) is replaced by the three things the
    mapper actually reads from it: cell connectivity, vertex positions, cell centres
  * per-cell 2D Kabsch, z-squish volumes and per-point barycentric volumes are computed
    as batched numpy/LAPACK calls instead of Python loops
  * pygcode is replaced by a regex reader producing the same floats (pygcode uses float(text))
  * find_closest_cell is only evaluated for points that have no containing cell
    (the notebook computes it for every point but only uses it for those)
  * the stateful smoothing/interpolation loop runs on plain Python floats
  * output lines are formatted with the same format specs and written in one go
"""
import math
import os
import re

import numpy as np

from .geometry import rotation_matrices, tangential_vectors  # noqa: F401 (re-exported: older imports)
from .timing import TIMER

# The notebook-exact baseline: every output fix off, the notebook's 1 mm retraction. map_gcode fills in from here
# whatever its caller leaves out. The user-facing defaults (fixes on, with what each flag does) are
# s4/profile.py MAP_DEFAULTS; the pipeline always passes a full map section from there.
MAPPING_DEFAULTS = dict(
    SEG_SIZE=0.6,
    MAX_ROTATION=30,
    MIN_ROTATION=-130,
    NOZZLE_OFFSET=42,
    ROTATION_AVERAGING_ALPHA=0.2,
    RETRACTION_LENGTH=1.0,
    ROTATION_MAX_DELTA=float(np.deg2rad(1)),
    MAX_EXTRUSION_MULTIPLIER=10,
    SPLIT_RETRACTIONS=False,
    SMOOTH_EXTRUSION_MULTIPLIER=False,
    EXTRUSION_MULTIPLIER_RANGE=None,
    SAFE_TRAVEL_TRANSITIONS=False,
)

_WORD = re.compile(r"([A-Za-z])\s*(-?(?:\d+\.?\d*|\.\d+))")


def _abs_det3(rows):
    """|det| of stacked 3x3 matrices, same LAPACK path as np.linalg.det on one matrix."""
    return np.abs(np.linalg.det(rows))


def tet_volumes(p1, p2, p3, p4):
    """Vectorised reference.tetrahedron_volume: |det(vstack([p2-p1, p3-p1, p4-p1]))| / 6."""
    return _abs_det3(np.stack([p2 - p1, p3 - p1, p4 - p1], axis=1)) / 6


def c_split_fractions(x0, y0, x1, y1, cstep, n_rot=1):
    """MAX_C_STEP (shared with s4/reference.py, so both cut at the same floats): where to cut a print move from
    (x0, y0) to (x1, y1) so that no piece sweeps more than cstep rad of C (C is the polar angle of the position).

    The cuts are where the rays at th0 + k * dth / n (k = 1 .. n - 1, n = int(|dth| / cstep) + 1) cross the move,
    so every piece sweeps the same angle, dth / n. (Cutting into equal lengths instead puts nearly all of the sweep
    of a move passing close to the axis into its middle piece.) A straight move that misses the axis turns
    monotonically by less than 180 deg, so each ray crosses it once. n_rot > 1 adds the cuts k / n_rot of the
    B-interpolation split (ROTATION_MAX_DELTA). Returns the increasing fractions of the move at which the pieces
    end (the last one 1.0), or None when the move sweeps at most cstep. A move along a line through the axis
    (C flips by 180 deg at r = 0, which no cut changes) is not cut."""
    th0 = math.atan2(y0, x0)
    dth = math.atan2(y1, x1) - th0
    dth = (dth + math.pi) % (2 * math.pi) - math.pi
    if abs(dth) <= cstep:
        return None
    dx, dy = x1 - x0, y1 - y0
    if x0 * dy - y0 * dx == 0:
        return None
    n = int(abs(dth) / cstep) + 1
    cuts = []
    for k in range(1, n):
        phi = th0 + dth * k / n
        ux, uy = math.cos(phi), math.sin(phi)
        den = ux * dy - uy * dx
        if den != 0:
            s = (uy * x0 - ux * y0) / den
            if 0.0 < s < 1.0:
                cuts.append(s)
    cuts.extend(k / n_rot for k in range(1, n_rot))
    out = []
    for s in sorted(cuts):
        if not out or s > out[-1]:
            out.append(s)
    out.append(1.0)
    return out


def split_pieces(fractions, extrusion, inv_time_feed):
    """Per piece of a move cut at `fractions` (c_split_fractions): (end fraction, E, inverse-time feed). A piece's E
    is its share of the move's E (the last piece takes the rest, so they sum to the move's E) and its time is its
    share of the move's time, so the move takes as long as before: inverse-time feed = the move's / share."""
    out = []
    prev = 0.0
    left = extrusion
    last = len(fractions) - 1
    for k, s in enumerate(fractions):
        share = s - prev
        prev = s
        if extrusion is None:
            e = None
        elif k == last:
            e = left
        else:
            e = extrusion * share
            left = left - e
        out.append((s, e, inv_time_feed / share if inv_time_feed is not None else None))
    return out


_POINTS_CACHE = {}  # (path, size, mtime, SEG_SIZE) -> read_gcode_points result (the mapper and the support check)


def read_gcode_points(path, SEG_SIZE):
    """Equivalent of reference.read_gcode_points. Returns dict of arrays/lists, plus "line": the index (in
    str.splitlines order) of the G-code line each point comes from. The result is cached per file and SEG_SIZE (the
    mapper and the support check read the same planar file); treat it as read-only."""
    st = os.stat(path)
    key = (os.path.abspath(path), st.st_size, st.st_mtime_ns, float(SEG_SIZE))
    hit = _POINTS_CACHE.get(key)
    if hit is None:
        _POINTS_CACHE.clear()
        hit = _POINTS_CACHE[key] = _read_gcode_points(path, SEG_SIZE)
    return hit


def _ffill(has, values, initial):
    """values[i] where has[i], else the last earlier value that had one (initial before the first)."""
    idx = np.where(has, np.arange(len(has)), -1)
    np.maximum.accumulate(idx, out=idx)
    out = np.where(idx >= 0, values[np.maximum(idx, 0)], initial)
    return out, idx


def _read_gcode_points(path, SEG_SIZE):
    """The parse is per line (regex, same words as before); the segmentation is vectorised with the same float
    operations as the per-line loop (_read_gcode_points_loop, kept as the reference): row-wise norm == 1-D norm,
    numpy floor division, and every quotient/product in the same order."""
    with open(path, "r") as fh:
        text = fh.read()
    line_no, mot = [], []
    vals = []  # per motion line: x, y, z, f, e (nan = absent)
    has = []
    nan = float("nan")
    for ln, raw in enumerate(text.splitlines()):
        code = raw.split(";", 1)[0]
        if "(" in code:
            code = re.sub(r"\([^)]*\)", "", code)
        words = _WORD.findall(code)
        if not words:
            continue
        motion = None
        x = y = z = f_val = e_val = None
        for letter, num in words:
            L = letter.upper()
            if L == "G":
                v = float(num)
                if v == 0:
                    motion = 0
                elif v == 1:
                    motion = 1
            elif L == "X":
                x = float(num)
            elif L == "Y":
                y = float(num)
            elif L == "Z":
                z = float(num)
            elif L == "F":
                f_val = float(num)
            elif L == "E":
                e_val = float(num)
        if motion is None:
            continue
        line_no.append(ln)
        mot.append(motion)
        vals.append((nan if x is None else x, nan if y is None else y, nan if z is None else z,
                     nan if f_val is None else f_val, nan if e_val is None else e_val))
        has.append((x is not None, y is not None, z is not None, f_val is not None, e_val is not None))
    m = len(mot)
    if m == 0:
        return {"position": np.zeros((0, 3)), "command": [], "extrusion": [], "inv_time_feed": [], "feed": [],
                "line": np.zeros(0, np.int64)}
    V = np.array(vals, dtype=np.float64)
    H = np.array(has, dtype=bool)
    pos = np.empty((m, 3))
    for k, init in enumerate((0.0, 0.0, 20.0)):
        pos[:, k] = _ffill(H[:, k], V[:, k], init)[0]
    prev = np.vstack([[0.0, 0.0, 20.0], pos[:-1]])
    delta = pos - prev
    distance = np.linalg.norm(delta, axis=1)
    # feed: sticky, starting at the int 5000 like the loop (kept as Python objects for the "feed" list)
    _, fidx = _ffill(H[:, 3], V[:, 3], 0.0)
    feed_obj = [5000] + V[:, 3].tolist()
    feed_at = [feed_obj[j + 1] for j in fidx.tolist()]
    feed_f = np.array(feed_at, dtype=np.float64)
    if (feed_f == 0).any():  # the loop (and the reference) divide by the feed: 1 / 0 raised there, here it gave F0
        ln = int(line_no[int(np.argmax(feed_f == 0))])
        raise ValueError(f"{path}: line {ln + 1}: a move at feed rate 0 (F0) has no finite time; give the move a feed rate")
    moving = distance > 0
    nseg = np.ones(m)
    nseg[moving] = -(-distance[moving] // SEG_SIZE)
    seg_distance = np.where(moving, distance / nseg, distance)
    t = (1 / feed_f) * seg_distance
    with np.errstate(divide="ignore"):
        inv = np.where(t == 0, np.nan, 1 / t)
    n = nseg.astype(np.int64)
    rep = np.repeat(np.arange(m), n)
    first = np.cumsum(n) - n
    k = (np.arange(len(rep)) - first[rep] + 1).astype(np.float64)
    mv = moving[rep]
    P = np.where(mv[:, None], prev[rep] + delta[rep] * k[:, None] / nseg[rep][:, None], pos[rep])
    e_has = H[:, 4]
    e_div = np.where(moving, V[:, 4] / nseg, V[:, 4])
    e_list = e_div.tolist()
    inv_list = inv.tolist()
    mot_s = ("G00", "G01")
    rep_l = rep.tolist()
    ehl = e_has.tolist()
    commands = [mot_s[mot[i]] for i in rep_l]
    extrusions = [e_list[i] if ehl[i] else None for i in rep_l]
    inv_ok = (~np.isnan(inv)).tolist()
    inv_feeds = [inv_list[i] if inv_ok[i] else None for i in rep_l]
    feeds = [feed_at[i] for i in rep_l]
    return {"position": P, "command": commands, "extrusion": extrusions, "inv_time_feed": inv_feeds,
            "feed": feeds, "line": np.asarray(line_no, np.int64)[rep]}


def _read_gcode_points_loop(path, SEG_SIZE):
    """Equivalent of reference.read_gcode_points. Returns dict of arrays/lists."""
    pos = np.array([0., 0., 20.])
    feed = 5000
    positions = []
    commands = []
    extrusions = []
    inv_feeds = []
    feeds = []
    with open(path, "r") as fh:
        text = fh.read()
    for ln, raw in enumerate(text.splitlines()):
        code = raw.split(";", 1)[0]
        if "(" in code:
            code = re.sub(r"\([^)]*\)", "", code)
        words = _WORD.findall(code)
        if not words:
            continue
        motion = None
        x = y = z = None
        f_val = None
        e_val = None
        for letter, num in words:
            L = letter.upper()
            if L == "G":
                v = float(num)
                if v == 0:
                    motion = "G00"
                elif v == 1:
                    motion = "G01"
            elif L == "X":
                x = float(num)
            elif L == "Y":
                y = float(num)
            elif L == "Z":
                z = float(num)
            elif L == "F":
                f_val = float(num)
            elif L == "E":
                e_val = float(num)
        if motion is None:
            continue
        prev_pos = pos.copy()
        if x is not None:
            pos[0] = x
        if y is not None:
            pos[1] = y
        if z is not None:
            pos[2] = z
        if f_val is not None:
            feed = f_val
        if feed == 0:
            raise ValueError(f"{path}: line {ln + 1}: a move at feed rate 0 (F0) has no finite time; give the move a feed rate")
        extrusion = e_val
        delta_pos = pos - prev_pos
        distance = np.linalg.norm(delta_pos)
        if distance > 0:
            num_segments = -(-distance // SEG_SIZE)
            seg_distance = distance / num_segments
            t = (1 / feed) * seg_distance
            inv = None if t == 0 else 1 / t
            n = int(num_segments)
            k = np.arange(1, n + 1, dtype=float)
            positions.append(prev_pos[None, :] + delta_pos[None, :] * k[:, None] / num_segments)
            commands.extend([motion] * n)
            e = extrusion / num_segments if extrusion is not None else None
            extrusions.extend([e] * n)
            inv_feeds.extend([inv] * n)
            feeds.extend([feed] * n)
        else:
            t = (1 / feed) * distance
            inv = None if t == 0 else 1 / t
            positions.append(pos.copy()[None, :])
            commands.append(motion)
            extrusions.append(extrusion)
            inv_feeds.append(inv)
            feeds.append(feed)
    P = np.concatenate(positions) if positions else np.zeros((0, 3))
    return {"position": P, "command": commands, "extrusion": extrusions, "inv_time_feed": inv_feeds, "feed": feeds}


def cell_rotations_and_squish(input_cells, input_points, input_centers, def_points, def_centers, MAX_ROTATION,
                              MIN_ROTATION):
    """Vectorised per-cell 2D Kabsch rotation, per-vertex averaged rotation, and z-squish scale."""
    n_cells = input_cells.shape[0]
    new_v = def_points[input_cells] - def_centers[:, None, :]
    old_v = input_points[input_cells] - input_centers[:, None, :]

    with np.errstate(invalid="ignore", divide="ignore"):
        px = input_centers[:, :2] / np.linalg.norm(input_centers[:, :2], axis=1)[:, None]
    plane_x = np.concatenate([px, np.zeros((n_cells, 1))], axis=1)

    def project(v):
        # reference: projected_x = np.sum(plane_x_axis * point, axis=1); projected_y likewise with [0,0,1]
        proj_x = np.sum(plane_x[:, None, :] * v, axis=2)
        proj_y = np.sum(np.array([0, 0, 1])[None, None, :] * v, axis=2)
        return np.stack([proj_x, proj_y], axis=2)  # (n, 4, 2)

    newp = project(new_v)
    oldp = project(old_v)
    cov = np.matmul(np.swapaxes(newp, 1, 2), oldp)  # (n, 2, 2)
    U, _, Vt = np.linalg.svd(cov)
    rot = np.matmul(U, Vt)
    r00 = np.clip(rot[:, 0, 0], -1, 1)
    rotation = -np.arccos(r00)
    rotation = np.where(rot[:, 1, 0] < 0, -rotation, rotation)
    rotation = np.maximum(np.minimum(rotation, np.deg2rad(MAX_ROTATION)), np.deg2rad(MIN_ROTATION))

    # vertex_rotations[v] += rotation / num_cells_per_vertex[v], accumulated in cell order
    counts = np.bincount(input_cells.ravel(), minlength=len(input_points)).astype(float)
    contrib = rotation[:, None] / counts[input_cells]
    vertex_rotations = np.zeros(len(def_points))
    np.add.at(vertex_rotations, input_cells.ravel(), contrib.ravel())

    unwarped = input_points[input_cells]
    warped = def_points[input_cells]
    vol0 = tet_volumes(*[unwarped[:, i] for i in range(4)])
    vold = tet_volumes(*[warped[:, i] for i in range(4)])
    z_squish = vol0 / vold
    return rotation, vertex_rotations, z_squish, vertex_volume_ratio(input_cells, vol0, vold, len(def_points))


def vertex_volume_ratio(cells, vol0, vold, n_points):
    """Per-vertex volume ratio: undeformed / deformed volume of all tets touching the vertex."""
    v0 = np.zeros(n_points)
    vd = np.zeros(n_points)
    np.add.at(v0, cells.ravel(), np.repeat(vol0, 4))
    np.add.at(vd, cells.ravel(), np.repeat(vold, 4))
    return v0 / vd


def map_gcode(input_cells, input_points, input_centers, def_points, def_centers, find_cells, gcode_path, mp=None):
    """
    find_cells(points) -> containing cell id per point (closest cell where none contains it), implemented by the
    caller with VTK (meshio_s4.make_find_cells). Returns (new_points dict-of-lists, stats).
    """
    mp = {**MAPPING_DEFAULTS, **(mp or {})}

    with TIMER("per-cell rotations + z squish (vectorised)"):
        _, vertex_rotations, z_squish_scales, vertex_ratio = cell_rotations_and_squish(
            input_cells, input_points, input_centers, def_points, def_centers, mp["MAX_ROTATION"], mp["MIN_ROTATION"])
        vertex_transformations = def_points - input_points

    with TIMER("read + segment gcode (regex)"):
        g = read_gcode_points(gcode_path, mp["SEG_SIZE"])
    P = g["position"]

    with TIMER("find containing/closest cells (VTK)"):
        containing = find_cells(P)

    with TIMER("barycentric (vectorised)"):
        idx = input_cells[containing]            # (N,4) vertex ids of containing cell
        cv = def_points[idx]                      # (N,4,3)
        a, b, c, d = cv[:, 0], cv[:, 1], cv[:, 2], cv[:, 3]
        total = tet_volumes(a, b, c, d)
        if np.any(total == 0):
            raise ValueError("The points do not form a valid tetrahedron (zero volume).")
        vol_a = tet_volumes(P, b, c, d)
        vol_b = tet_volumes(P, a, c, d)
        vol_c = tet_volumes(P, a, b, d)
        vol_d = tet_volumes(P, a, b, c)
        bary = np.stack([vol_a / total, vol_b / total, vol_c / total, vol_d / total], axis=1)
        bary_ok = ~(np.sum(bary, axis=1) > 1.01)
        transformation = np.sum(vertex_transformations[idx] * bary[:, :, None], axis=1)
        new_pos_all = P - transformation
        rot_all = np.sum(vertex_rotations[idx] * bary, axis=1)
        if mp["SMOOTH_EXTRUSION_MULTIPLIER"]:
            squish_all = np.sum(vertex_ratio[idx] * bary, axis=1)
        else:
            squish_all = z_squish_scales[containing]
    smooth = map_smooth_of(mp)
    if smooth is not None:
        with TIMER("smooth mapping (MAP_SMOOTH)"):
            new_pos_all = smooth_mapping(P, new_pos_all, bary_ok, g, mp["RETRACTION_LENGTH"], *smooth)
    # for the support check (s4/support_check.py check(mapped=...)): the positions the G-code is written from
    LAST_MAP.clear()
    LAST_MAP.update(planar=gcode_path, seg_size=mp["SEG_SIZE"], g=g, real=new_pos_all, ok=bary_ok,
                    map_smooth=smooth, retraction=mp["RETRACTION_LENGTH"])

    with TIMER("sequential smoothing loop"):
        out = _sequential(new_pos_all, rot_all, bary_ok, squish_all, g, mp)
    stats = {
        "lost_vertices": out.pop("_lost"),
        "gcode_points": len(P),
        "no_containing_cell": 0,
        "bad_barycentric_sum": int((~bary_ok).sum()),
    }
    return out, stats


SMOOTH_STATS = {}  # diagnostics of the last smooth_mapping call
LAST_MAP = {}  # the last map_gcode's planar points, their real-space positions (after MAP_SMOOTH) and in-mesh flags


def map_smooth_of(mp):
    """(sigma, max_shift) of MAP_SMOOTH, or None when it is off: what smooth_mapping gets from the map settings."""
    if mp.get("MAP_SMOOTH", 0) and mp["MAP_SMOOTH"] > 0:
        return float(mp["MAP_SMOOTH"]), float(mp.get("MAP_SMOOTH_MAX", 0.2))
    return None


def smooth_mapping(P, mapped, ok, g, retraction, sigma, max_shift, window=12, chunk=100000):
    """MAP_SMOOTH: a smoother map for the printed paths (opt-in; s4/reference.py calls this same function).

    The deformation is linear per tet, so a straight planar line maps to a polyline that kinks wherever it crosses a
    tet face. Here the displacement D = P - mapped is re-estimated at every printed point by a local linear fit
    D ~ a + B (xy - xy_i) over its neighbours along the same printed run (same layer, no travel in between),
    Gaussian-weighted by in-plane distance (sigma mm). A linear fit reproduces an affine map exactly, so curves and
    corners of the planar path keep their shape; only the jumps of the gradient at tet faces are blended over about
    sigma. The change per point is capped at max_shift mm. Points outside the mesh (not ok) are neither used nor
    moved."""
    n = len(P)
    if n < 3:
        return mapped
    ext = np.array([np.nan if e is None else e for e in g["extrusion"]], dtype=float)
    printing = np.array([c == "G01" for c in g["command"]]) & (ext > 0) & (ext != retraction)
    member = printing | np.r_[printing[1:], False]  # printed points and the start point of each printed run
    run = np.cumsum(~printing)  # a run's start point (a travel end) shares the id of the moves printed from it
    D = P - mapped
    out = D.copy()
    use = member & ok
    offs = np.arange(-window, window + 1)
    xy = P[:, :2]
    inv2s2 = 0.5 / sigma ** 2
    for c0 in range(0, n, chunk):
        i = np.arange(c0, min(n, c0 + chunk))
        i = i[use[i]]
        if not len(i):
            continue
        j = i[:, None] + offs[None, :]
        inr = (j >= 0) & (j < n)
        j = np.clip(j, 0, n - 1)
        valid = inr & use[j] & (run[j] == run[i][:, None]) & (P[j, 2] == P[i, 2][:, None])
        dx = xy[j] - xy[i][:, None, :]
        w = np.exp(-np.einsum("ijk,ijk->ij", dx, dx) * inv2s2) * valid
        X = np.concatenate([np.ones(dx.shape[:2] + (1,)), dx], axis=2)  # (m, k, 3)
        A = np.einsum("ij,ijk,ijl->ikl", w, X, X)
        lam = 1e-3 * sigma ** 2 * w.sum(1) + 1e-12
        A[:, 1, 1] += lam
        A[:, 2, 2] += lam
        b = np.einsum("ij,ijk,ijl->ikl", w, X, D[j])
        a = np.linalg.solve(A, b)[:, 0, :]
        out[i] = a
    shift = out - D
    s = np.linalg.norm(shift, axis=1)
    big = s > max_shift
    shift[big] *= (max_shift / s[big])[:, None]
    SMOOTH_STATS.clear()
    SMOOTH_STATS.update(points=int(use.sum()), mean_shift_mm=round(float(s[use].mean()), 4) if use.any() else 0.0,
                        p99_shift_mm=round(float(np.percentile(s[use], 99)), 4) if use.any() else 0.0,
                        capped=int(big.sum()))
    return mapped - shift


def _sequential(new_pos_all, rot_all, bary_ok, squish_all, g, mp):
    """Walk the segmented planar moves in order and emit the 4-axis moves (notebook cell 17's loop).

    Inputs per planar point: new_pos_all = position mapped back to the part, rot_all = B tilt (rad), bary_ok =
    whether it lies in a tet, squish_all = extrusion multiplier; g = the planar moves (commands, E, feeds).
    Per point: points outside the part drop print moves and turn travels into a lift to the highest printed
    point; B is smoothed (ROTATION_AVERAGING_ALPHA) and steps over ROTATION_MAX_DELTA are split; E is scaled.
    The o_* lists are the output moves, one entry each: position, rotation, command, extrusion, inverse-time feed,
    travelling (z-hop) flag, planar feed, e_only (zero-motion retract/unretract line).

    Flags (all False = notebook-exact; S4_PIPELINE.md 'Differences from the notebook'): SPLIT (SPLIT_RETRACTIONS)
    retracts/unretracts in place as e_only moves; SAFE (SAFE_TRAVEL_TRANSITIONS) lowers to the true re-entry
    height after a travel that left the part, keeps each split step's own command and clamps z < 0 to the bed;
    MRANGE (EXTRUSION_MULTIPLIER_RANGE) clamps the multiplier. Returns the o_* lists as a dict, plus '_lost'.
    """
    ALPHA = mp["ROTATION_AVERAGING_ALPHA"]
    RET = mp["RETRACTION_LENGTH"]
    MAXD = mp["ROTATION_MAX_DELTA"]
    CSTEP = math.radians(mp.get("MAX_C_STEP", 0) or 0)  # 0 = off
    MAXE = mp["MAX_EXTRUSION_MULTIPLIER"]
    SPLIT = mp["SPLIT_RETRACTIONS"]
    MRANGE = mp["EXTRUSION_MULTIPLIER_RANGE"]
    SAFE = mp["SAFE_TRAVEL_TRANSITIONS"]
    lim45 = float(np.deg2rad(45))

    commands = g["command"]
    extrusions = g["extrusion"]
    inv_feeds = g["inv_time_feed"]
    feeds = g["feed"]
    newp = new_pos_all.tolist()
    rots = rot_all.tolist()
    ok = bary_ok.tolist()
    squish = squish_all.tolist()

    o_pos, o_rot, o_cmd, o_ext, o_inv, o_trav, o_feed, o_eonly = [], [], [], [], [], [], [], []
    o_src, o_kind = [], []  # per output move: planar point index, and kind (diagnostics only, never written)
    cur = [0, 0]

    def emit(pos, rot, cmd, ext, inv, trav, feed, e_only=False):
        o_src.append(cur[0])
        o_kind.append(3 if e_only else cur[1])
        o_pos.append(pos)
        o_rot.append(rot)
        o_cmd.append(cmd)
        o_ext.append(ext)
        o_inv.append(inv)
        o_trav.append(trav)
        o_feed.append(feed)
        o_eonly.append(e_only)
    prev_new_position = None
    travelling_over_air = False
    travelling = False
    prev_rotation = 0
    prev_travelling = False
    prev_command = "G00"
    highest_printed_point = 0
    lost = 0

    for i in range(len(commands)):
        cur[0], cur[1] = i, 0
        command = commands[i]
        inv_time_feed = inv_feeds[i]
        extrusion = extrusions[i]
        dont_smooth_rotation = False
        reentry = None
        if ok[i]:
            new_position = list(newp[i])
            rotation = rots[i]
            if SAFE and new_position[2] < 0:
                # a point just outside the part can map a little below the bed; write_gcode drops z < 0, which
                # lost the travel down to the re-entry point and left a print move extruding straight down
                new_position[2] = 0.0
            if travelling_over_air:
                if SAFE:
                    reentry = list(new_position)
                    new_position[2] = max(highest_printed_point, new_position[2])
                else:
                    new_position[2] = highest_printed_point
                rotation = max(min(rotation, lim45), -lim45)
                dont_smooth_rotation = True
            travelling_over_air = False
        else:
            if command == "G01":
                lost += 1
                continue
            elif command == "G00" and not travelling_over_air and prev_new_position is not None:
                new_position = [prev_new_position[0], prev_new_position[1], highest_printed_point]
                rotation = max(min(prev_rotation, lim45), -lim45)
                dont_smooth_rotation = True
                travelling_over_air = True
            else:
                continue

        move_command = command
        if reentry is not None and command == "G01" and extrusion is not None and extrusion != RET and extrusion != -RET:
            # a print move whose start was off the part: travel over its end point instead of extruding in the air
            move_command = "G00"; extrusion = None; lost += 1
        extrusion_multiplier = 1
        if extrusion is not None and extrusion != RET and extrusion != -RET:
            extrusion_multiplier = extrusion_multiplier * squish[i]
            if MRANGE is not None:
                extrusion_multiplier = min(max(extrusion_multiplier, MRANGE[0]), MRANGE[1])
            extrusion = extrusion * min(extrusion_multiplier, MAXE)
        elif extrusion == -RET:
            travelling = True
        elif extrusion == RET:
            travelling = False
        if prev_rotation is not None and not dont_smooth_rotation:
            rotation = ALPHA * rotation + (1 - ALPHA) * prev_rotation

        # retraction event: optionally move the E value onto its own zero-motion line
        split = SPLIT and extrusion is not None and (extrusion == RET or extrusion == -RET)
        motion_extrusion = None if split else extrusion
        if split and extrusion < 0:  # retract in place first, then lift
            if prev_new_position is not None:
                emit(tuple(prev_new_position), prev_rotation, "G01", extrusion, None, prev_travelling, feeds[i], True)
            else:
                emit(tuple(new_position), rotation, "G01", extrusion, None, prev_travelling, feeds[i], True)

        # MAX_C_STEP: the machine moves C and X linearly, so near the turntable axis a move that sweeps a large C
        # angle is drawn as a curve; print moves (not the unretract) are cut so every piece sweeps at most CSTEP
        rot_split = prev_new_position is not None and abs(rotation - prev_rotation) > MAXD
        fr = None
        if (CSTEP and prev_new_position is not None and motion_extrusion is not None and motion_extrusion > 0
                and motion_extrusion != RET and motion_extrusion != -RET):
            fr = c_split_fractions(prev_new_position[0], prev_new_position[1], new_position[0], new_position[1],
                                   CSTEP, int(abs(rotation - prev_rotation) / MAXD) + 1 if rot_split else 1)
        if rot_split:
            delta_rotation = rotation - prev_rotation
            n = int(abs(delta_rotation) / MAXD) + 1
            dx = new_position[0] - prev_new_position[0]
            dy = new_position[1] - prev_new_position[1]
            dz = new_position[2] - prev_new_position[2]
            cur[1] = 1  # rotation-split steps
            if fr is None:
                e_i = motion_extrusion / n if motion_extrusion is not None else None
                f_i = inv_time_feed * n if inv_time_feed is not None else None
                for k in range(n):
                    s = (k + 1) / n
                    emit((prev_new_position[0] + dx * s, prev_new_position[1] + dy * s, prev_new_position[2] + dz * s),
                         prev_rotation + delta_rotation * s, move_command if SAFE else prev_command, e_i, f_i,
                         prev_travelling, feeds[i])
            else:  # the B steps' cuts and the C steps' cuts together
                for s, e_i, f_i in split_pieces(fr, motion_extrusion, inv_time_feed):
                    emit((prev_new_position[0] + dx * s, prev_new_position[1] + dy * s, prev_new_position[2] + dz * s),
                         prev_rotation + delta_rotation * s, move_command if SAFE else prev_command, e_i, f_i,
                         prev_travelling, feeds[i])
        elif fr is not None:
            cur[1] = 4  # C-split steps
            dx = new_position[0] - prev_new_position[0]
            dy = new_position[1] - prev_new_position[1]
            dz = new_position[2] - prev_new_position[2]
            dr = rotation - prev_rotation
            for s, e_i, f_i in split_pieces(fr, motion_extrusion, inv_time_feed):
                emit((prev_new_position[0] + dx * s, prev_new_position[1] + dy * s, prev_new_position[2] + dz * s),
                     prev_rotation + dr * s, move_command, e_i, f_i, travelling, feeds[i])
        else:
            emit(tuple(new_position), rotation, move_command, motion_extrusion, inv_time_feed, travelling, feeds[i])

        cur[1] = 2
        if reentry is not None and reentry[2] != new_position[2]:  # lower to the true height, as a travel
            emit(tuple(reentry), rotation, "G00", None, None, travelling, feeds[i])
            new_position = reentry

        if split and extrusion > 0:  # plunge first, then unretract in place
            if o_trav[-1] != travelling:  # rotation-interpolated motion stayed hopped: plunge without E
                emit(tuple(new_position), rotation, command, None, None, travelling, feeds[i])
            emit(tuple(new_position), rotation, "G01", extrusion, None, travelling, feeds[i], True)

        prev_rotation = rotation
        prev_new_position = list(new_position)
        prev_travelling = travelling
        prev_command = move_command if SAFE else command

        if (command == "G01" and extrusion is not None and extrusion > 0
                and (highest_printed_point != 0 or new_position[2] < 1)):
            highest_printed_point = max(highest_printed_point, new_position[2])

    return {"position": o_pos, "rotation": o_rot, "command": o_cmd, "extrusion": o_ext,
            "inv_time_feed": o_inv, "travelling": o_trav, "feed": o_feed, "e_only": o_eonly, "_lost": lost,
            "src": o_src, "kind": o_kind}


def write_gcode(pts, out_path, NOZZLE_OFFSET=42):
    """Equivalent of reference.write_gcode (cell 18)."""
    P = np.asarray(pts["position"], dtype=float).reshape(-1, 3)
    rot = np.asarray(pts["rotation"], dtype=float)
    keep = ~(np.all(np.isnan(P), axis=1)) & ~(P[:, 2] < 0)
    idx = np.nonzero(keep)[0]
    P = P[idx]
    rot = rot[idx]
    trav = np.asarray(pts["travelling"], dtype=bool)[idx]
    z_hop = trav.astype(float)

    r = np.linalg.norm(P[:, :2], axis=1)
    theta = np.arctan2(P[:, 1], P[:, 0])
    z = P[:, 2].copy()
    r += -np.sin(rot) * (NOZZLE_OFFSET + z_hop)
    z += (np.cos(rot) - 1) * (NOZZLE_OFFSET + z_hop) + z_hop

    prev = np.concatenate([[0.0], theta[:-1]])
    delta = theta - prev
    delta = np.where(delta > np.pi, delta - 2 * np.pi, delta)
    delta = np.where(delta < -np.pi, delta + 2 * np.pi, delta)
    theta_accum = np.cumsum(delta)
    C = np.rad2deg(theta_accum)
    B = np.rad2deg(rot)

    cmds = pts["command"]
    exts = pts["extrusion"]
    invs = pts["inv_time_feed"]
    e_only = pts.get("e_only") or [False] * len(cmds)
    feeds = pts.get("feed")
    lines = [
        "G94 ; mm/min feed  \n",
        "G28 ; home \n",
        "M83 ; relative extrusion \n",
        "G1 E10 ; prime extruder \n",
        "G94 ; mm/min feed \n",
        "G90 ; absolute positioning \n",
        "G0 C0 X0 Z20 B0 ; go to start \n",
        "G93 ; inverse time feed \n",
    ]
    Cl, rl, zl, Bl = C.tolist(), r.tolist(), z.tolist(), B.tolist()
    for j, i in enumerate(idx.tolist()):
        s = f"{cmds[i]} C{Cl[j]:.5f} X{rl[j]:.5f} Z{zl[j]:.5f} B{Bl[j]:.5f}"
        e = exts[i]
        if e is not None:
            s += f" E{e:.4f}"
        f = invs[i]
        if e_only[i]:  # zero-motion retract/unretract at the planar G-code's feed (mm/min)
            lines.append("G94\n")
            lines.append(s + f" F{feeds[i]:g}\n")
            lines.append("G93\n")
        elif f is not None:
            lines.append(s + f" F{f:.4f}\n")
        else:
            lines.append("G94\n")
            lines.append(s + " F20000\n")
            lines.append("G93\n")
    with open(out_path, "w") as fh:
        fh.write("".join(lines))
    return idx  # the output moves written, one G00/G01 line each, in order
