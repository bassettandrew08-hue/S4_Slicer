"""
Quality checks of a final 4-axis G-code.

vertical_extrusion finds poles and extrusion along the nozzle's own axis (the pipeline prints it as [quality]).

path_roughness, for comparing runs (tools/model_suite.py), is R-Theta Sim's "Path roughness" view in Python: the
extruding moves are split into bead segments the way the sim does (by C every 3 deg and B every 6 deg, at most 32),
and each segment gets the angle it turns from the segment it continues. A turn only counts when a neighbouring
segment also turns sharply: a real corner of the part is one sharp turn between straight runs, while zig-zag skin
and jittery mapping turn sharply again and again.
"""
import re

import numpy as np

from .geometry import nozzle_tip

_WORD = re.compile(r"([CXZBEF])(-?\d+(?:\.\d*)?|-?\.\d+)")


def bead_segments(gcode_path, nozzle_offset=42.0):
    """Bead segments of the extruding moves: (start xyz, end xyz) in the bed frame, z up. Same split as the sim."""
    C = X = Z = B = 0.0
    rows = []  # per extruding move: C0 X0 Z0 B0 C1 X1 Z1 B1
    with open(gcode_path, errors="replace") as fh:
        for line in fh:
            if not line.startswith(("G0", "G1")):
                continue
            w = dict(_WORD.findall(line.split(";", 1)[0]))
            c, x, z, b = (float(w.get(k, v)) for k, v in (("C", C), ("X", X), ("Z", Z), ("B", B)))
            if line.startswith(("G1 ", "G01")) and float(w.get("E", 0)) > 1e-6 and (c, x, z, b) != (C, X, Z, B):
                rows.append((C, X, Z, B, c, x, z, b))
            C, X, Z, B = c, x, z, b
    if not rows:
        return np.zeros((0, 3)), np.zeros((0, 3))
    a = np.array(rows)
    k = np.clip(np.maximum(np.ceil(np.abs(a[:, 4] - a[:, 0]) / 3.0), np.ceil(np.abs(a[:, 7] - a[:, 3]) / 6.0)), 1, 32)
    k = k.astype(int)
    mv = np.repeat(np.arange(len(a)), k)
    j = np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)
    f0 = j / k[mv]
    f1 = (j + 1) / k[mv]

    def tip(f):
        s = a[mv, :4] + (a[mv, 4:] - a[mv, :4]) * f[:, None]
        x, y, h = nozzle_tip(s[:, 0], s[:, 1], s[:, 2], s[:, 3], nozzle_offset)
        return np.c_[x, y, h]
    return tip(f0), tip(f1)


def path_roughness(gcode_path, nozzle_offset=42.0, threshold=30.0):
    """Per-segment roughness (deg) and the count above `threshold`. Returns dict with arrays for region queries."""
    p0, p1 = bead_segments(gcode_path, nozzle_offset)
    n = len(p0)
    raw = np.zeros(n)
    if n > 1:
        d = p1 - p0
        L = np.linalg.norm(d, axis=1)
        cont = (np.linalg.norm(p0[1:] - p1[:-1], axis=1) <= 0.05) & (L[1:] > 1e-6) & (L[:-1] > 1e-6)
        cosang = np.einsum("ij,ij->i", d[1:], d[:-1]) / np.maximum(L[1:] * L[:-1], 1e-12)
        raw[1:] = np.where(cont, np.degrees(np.arccos(np.clip(cosang, -1, 1))), 0.0)
    nb = np.maximum(np.r_[0.0, raw[:-1]], np.r_[raw[1:], 0.0])
    rough = np.minimum(raw, nb)
    return {"segments": n, "rough": rough, "mid": (p0 + p1) / 2 if n else p0, "count": int((rough > threshold).sum())}


def vertical_extrusion(gcode_path, nozzle_offset=42.0, min_dz=1.0, pole_dz=2.0, retraction=1.0):
    """Suspicious extruding moves (longer than min_dz) in the final 4-axis G-code. Returns (poles, along_axis):
    poles start right after a travel and drop straight down more than pole_dz (extrusion dragged down from a lifted
    travel point: these print as free-standing sticks); along_axis ones run mostly along the nozzle's own axis
    (pushing into or pulling out of the bead). A vertical move with the nozzle tilted sideways is normal S4 printing
    and isn't flagged."""
    word = re.compile(r"([CXZBE])(-?\d+(?:\.\d*)?|-?\.\d+)")
    rows = []
    with open(gcode_path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.startswith(("G0", "G1")):
                continue
            w = dict(word.findall(line))
            if "X" in w:
                rows.append((n, float(w["C"]), float(w["X"]), float(w["Z"]), float(w["B"]), float(w.get("E", "nan"))))
    if len(rows) < 2:
        return [], []
    a = np.array(rows)
    b = np.radians(a[:, 4])
    th = np.radians(a[:, 1])
    x, y, z = nozzle_tip(a[:, 1], a[:, 2], a[:, 3], a[:, 4], nozzle_offset)
    P = np.c_[x, y, z]
    d = np.diff(P, axis=0)
    dz = np.abs(d[:, 2])
    dxy = np.hypot(d[:, 0], d[:, 1])
    E = a[:, 5]
    printing = (np.nan_to_num(E) > 0) & ~np.isclose(E, retraction)  # an unretract is not printing
    n = np.linalg.norm(d, axis=1)
    bb, tt = b[1:], th[1:]  # nozzle axis at the segment end: radial -sin B, vertical cos B
    axis = np.c_[-np.sin(bb) * np.cos(tt), -np.sin(bb) * np.sin(tt), np.cos(bb)]
    along = np.abs(np.sum(d * axis, axis=1)) > 0.894 * n  # within ~27 deg of the nozzle axis
    vertical = (dz > min_dz) & (dxy < 0.5 * dz)
    poles, steep = [], []
    for k in np.nonzero(printing[1:] & (n > min_dz) & (vertical | along))[0]:
        item = (int(a[k + 1, 0]), float(P[k, 2]), float(P[k + 1, 2]), [round(float(v), 1) for v in P[k + 1, :2]])
        if vertical[k] and not printing[k] and dz[k] > pole_dz:
            poles.append(item)
        elif along[k]:
            steep.append(item)
    return poles, steep


def format_vertical(result, line_offset=0):
    poles, steep = result
    s = f"[quality] poles (extruding >2 mm straight down from a travel): {len(poles) or 'none'}"
    for line, z0, z1, xy in poles[:3]:
        s += f"\n[quality]   G-code line {line + line_offset}: z {z0:.1f} -> {z1:.1f} at {xy}"
    if steep:
        s += ("\n[quality] extruding along the nozzle axis (>1 mm; pushing into or pulling out of the bead): "
              f"{len(steep)}")
    return s
