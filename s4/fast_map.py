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
import re

import numpy as np

from .timing import TIMER

MAPPING_DEFAULTS = dict(
    SEG_SIZE=0.6,
    MAX_ROTATION=30,
    MIN_ROTATION=-130,
    NOZZLE_OFFSET=42,
    ROTATION_AVERAGING_ALPHA=0.2,
    RETRACTION_LENGTH=1.0,
    ROTATION_MAX_DELTA=float(np.deg2rad(1)),
    MAX_EXTRUSION_MULTIPLIER=10,
    # False = notebook-exact. True = retract/unretract in place (E-only line at Cura's feed) instead of
    # extruding +/-RETRACTION_LENGTH during the 1 mm travel lift/plunge (which leaves filament "sticks").
    SPLIT_RETRACTIONS=False,
    # False = notebook-exact: extrusion is scaled by each tet's volume ratio (constant per tet, so it jumps at
    # tet boundaries). True = volume-weighted per-vertex ratio, interpolated barycentrically like position.
    SMOOTH_EXTRUSION_MULTIPLIER=False,
    EXTRUSION_MULTIPLIER_RANGE=None,  # e.g. (0.5, 2.0): clamp the volume-ratio multiplier
)

_WORD = re.compile(r"([A-Za-z])\s*(-?(?:\d+\.?\d*|\.\d+))")


def tangential_vectors(cell_centers):
    cxy = cell_centers[:, :2]
    c3 = np.hstack([cxy, np.zeros((cxy.shape[0], 1))])
    t = np.cross(np.array([0, 0, 1]), c3)
    with np.errstate(invalid="ignore"):
        t /= np.linalg.norm(t, axis=1)[:, None]
    t[np.isnan(t).any(axis=1)] = [1, 0, 0]
    return t


def rotation_matrices(cell_centers, rotation_field):
    from scipy.spatial.transform import Rotation as R
    return R.from_rotvec(rotation_field[:, None] * tangential_vectors(cell_centers)).as_matrix()


def _abs_det3(rows):
    """|det| of stacked 3x3 matrices, same LAPACK path as np.linalg.det on one matrix."""
    return np.abs(np.linalg.det(rows))


def tet_volumes(p1, p2, p3, p4):
    """Vectorised reference.tetrahedron_volume: |det(vstack([p2-p1, p3-p1, p4-p1]))| / 6."""
    return _abs_det3(np.stack([p2 - p1, p3 - p1, p4 - p1], axis=1)) / 6


def read_gcode_points(path, SEG_SIZE):
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
    for raw in text.splitlines():
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


def cell_rotations_and_squish(input_cells, input_points, input_centers, def_points, def_centers, MAX_ROTATION, MIN_ROTATION):
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
    find_cells(points) -> (containing, closest_or_None_fn) implemented by the caller with VTK.
    Returns (new_points dict-of-lists, stats).
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

    with TIMER("sequential smoothing loop"):
        out = _sequential(P, new_pos_all, rot_all, bary_ok, squish_all, g, mp)
    stats = {
        "lost_vertices": out.pop("_lost"),
        "gcode_points": len(P),
        "no_containing_cell": 0,
        "bad_barycentric_sum": int((~bary_ok).sum()),
    }
    return out, stats


def _sequential(P, new_pos_all, rot_all, bary_ok, squish_all, g, mp):
    ALPHA = mp["ROTATION_AVERAGING_ALPHA"]
    RET = mp["RETRACTION_LENGTH"]
    MAXD = mp["ROTATION_MAX_DELTA"]
    MAXE = mp["MAX_EXTRUSION_MULTIPLIER"]
    SPLIT = mp["SPLIT_RETRACTIONS"]
    MRANGE = mp["EXTRUSION_MULTIPLIER_RANGE"]
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

    def emit(pos, rot, cmd, ext, inv, trav, feed, e_only=False):
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
        command = commands[i]
        inv_time_feed = inv_feeds[i]
        extrusion = extrusions[i]
        dont_smooth_rotation = False
        if ok[i]:
            new_position = list(newp[i])
            rotation = rots[i]
            if travelling_over_air:
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

        if prev_new_position is not None and abs(rotation - prev_rotation) > MAXD:
            delta_rotation = rotation - prev_rotation
            n = int(abs(delta_rotation) / MAXD) + 1
            dx = new_position[0] - prev_new_position[0]
            dy = new_position[1] - prev_new_position[1]
            dz = new_position[2] - prev_new_position[2]
            e_i = motion_extrusion / n if motion_extrusion is not None else None
            f_i = inv_time_feed * n if inv_time_feed is not None else None
            for k in range(n):
                s = (k + 1) / n
                emit((prev_new_position[0] + dx * s, prev_new_position[1] + dy * s, prev_new_position[2] + dz * s),
                     prev_rotation + delta_rotation * s, prev_command, e_i, f_i, prev_travelling, feeds[i])
        else:
            emit(tuple(new_position), rotation, command, motion_extrusion, inv_time_feed, travelling, feeds[i])

        if split and extrusion > 0:  # plunge first, then unretract in place
            if o_trav[-1] != travelling:  # rotation-interpolated motion stayed hopped: plunge without E
                emit(tuple(new_position), rotation, command, None, None, travelling, feeds[i])
            emit(tuple(new_position), rotation, "G01", extrusion, None, travelling, feeds[i], True)

        prev_rotation = rotation
        prev_new_position = list(new_position)
        prev_travelling = travelling
        prev_command = command

        if command == "G01" and extrusion is not None and extrusion > 0 and (highest_printed_point != 0 or new_position[2] < 1):
            highest_printed_point = max(highest_printed_point, new_position[2])

    return {"position": o_pos, "rotation": o_rot, "command": o_cmd, "extrusion": o_ext,
            "inv_time_feed": o_inv, "travelling": o_trav, "feed": o_feed, "e_only": o_eonly, "_lost": lost}


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
