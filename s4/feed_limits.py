"""
Axis speed limits for the final 4-axis G-code (post-processing, shared by the fast and reference pipelines).

The notebook gives every move the time of its planar move (planar length / Cura's feed, as G93 inverse time). That
keeps the nozzle tip at Cura's speed, but ignores what the machine axes have to do to get it there: a 1 deg B change
swings X/Z by ~0.7 mm through the nozzle offset, a small move near the centre can need a big C rotation, and a travel
lifted over the part is much longer in real space than in the planar file. Those moves asked for thousands of deg/s
on C and B.

Each move's time becomes the longest of: the planned time, the tip distance at the move's feed, and each axis' travel
at its speed limit. Moves are only ever slowed down, never sped up. The G94 "F20000" moves (travel down to a re-entry
point) become G93 moves too, so every motion has a definite time.
"""
import re

import numpy as np

_WORD = re.compile(r"([A-Z])(-?\d+(?:\.\d*)?|-?\.\d+)")


def _fmt_f(f):
    return f"{f:.4f}" if f >= 1 else f"{f:.6g}"


def apply(path, max_c, max_b, max_x, max_z, nozzle_offset=42.0):
    """Rewrite the feeds in path in place. Speeds in deg/s (C, B) and mm/s (X, Z). Returns a summary dict."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    lim = {"C": max_c * 60.0, "B": max_b * 60.0, "X": max_x * 60.0, "Z": max_z * 60.0}  # per minute
    pos = {"C": 0.0, "X": 0.0, "Z": 20.0, "B": 0.0}
    out = []
    mode = 94
    started = False
    total_min = 0.0
    slowed = 0
    moves = 0
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        code = line.split(";", 1)[0].strip()
        if code == "G93":
            mode = 93
        elif code == "G94":
            mode = 94
            # a G94 / move F20000 / G93 triple (the writer's feedless move): make it a G93 move if it moves
            if i + 2 < n and lines[i + 2].strip() == "G93" and started and lines[i + 1].startswith(("G00 C", "G01 C")) \
                    and lines[i + 1].rstrip().endswith(" F20000"):
                w = dict(_WORD.findall(lines[i + 1]))
                new = {k: float(w[k]) for k in pos}
                if any(new[k] != pos[k] for k in pos):
                    t = _move_time(pos, new, None, 20000.0, lim, nozzle_offset)
                    out.append(lines[i + 1].rstrip()[:-len(" F20000")] + " F" + _fmt_f(1.0 / t))
                    pos = new; total_min += t; moves += 1; slowed += 1
                    mode = 93  # the triple's closing G93 is consumed here
                    i += 3
                    continue
        elif code.startswith(("G00 C", "G01 C", "G0 C", "G1 C")):
            w = dict(_WORD.findall(code))
            new = {k: float(w.get(k, pos[k])) for k in pos}
            if not started:  # "G0 C0 X0 Z20 B0 ; go to start"
                started = True
                pos = new
                out.append(line); i += 1
                continue
            if mode == 93 and "F" in w and any(new[k] != pos[k] for k in pos):
                t0 = 1.0 / float(w["F"])
                t = _move_time(pos, new, t0, None, lim, nozzle_offset)
                moves += 1; total_min += t
                if t > t0 * (1 + 1e-9):
                    slowed += 1
                    line = re.sub(r" F-?[\d.eE+-]+", " F" + _fmt_f(1.0 / t), line.rstrip())
            pos = new
        out.append(line)
        i += 1
    with open(path, "w") as fh:
        fh.write("\n".join(out))
    return {"moves_slowed_pct": round(100.0 * slowed / max(moves, 1), 1),
            "estimated_print_time_min": round(total_min, 1)}


def _tip(p, L):
    b = np.radians(p["B"]); th = np.radians(p["C"])
    r = p["X"] + np.sin(b) * L; z = p["Z"] - (np.cos(b) - 1) * L
    return np.array([r * np.cos(th), r * np.sin(th), z])


def _move_time(a, b, t0, feed, lim, L):
    """Minutes for the move a -> b."""
    t = t0 or 0.0
    if feed:
        t = max(t, float(np.linalg.norm(_tip(b, L) - _tip(a, L))) / feed)
    for k, v in lim.items():
        if v > 0:
            t = max(t, abs(b[k] - a[k]) / v)
    return max(t, 1e-9)
