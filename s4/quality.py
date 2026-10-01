"""
Surface-quality measures of a final 4-axis G-code, for comparing runs (tools/model_suite.py).

path_roughness is R-Theta Sim's "Path roughness" view in Python: the extruding moves are split into bead segments the
way the sim does (by C every 3 deg and B every 6 deg, at most 32), and each segment gets the angle it turns from the
segment it continues. A turn only counts when a neighbouring segment also turns sharply: a real corner of the part is
one sharp turn between straight runs, while zig-zag skin and jittery mapping turn sharply again and again.
"""
import re

import numpy as np

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
    f0 = j / k[mv]; f1 = (j + 1) / k[mv]

    def tip(f):
        s = a[mv, :4] + (a[mv, 4:] - a[mv, :4]) * f[:, None]
        bb = np.radians(s[:, 3]); th = np.radians(s[:, 0])
        r = s[:, 1] + np.sin(bb) * nozzle_offset
        h = s[:, 2] - (np.cos(bb) - 1) * nozzle_offset
        return np.c_[r * np.cos(th), r * np.sin(th), h]
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
