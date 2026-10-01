"""
Axis speed limits for the final 4-axis G-code (post-processing, shared by the fast and reference pipelines).

The notebook gives every move the time of its planar move (planar length / Cura's feed, as G93 inverse time). That
keeps the nozzle tip at Cura's speed, but ignores what the machine axes have to do to get it there: a 1 deg B change
swings X/Z by ~0.7 mm through the nozzle offset, a small move near the centre can need a big C rotation, and a travel
lifted over the part is much longer in real space than in the planar file. Those moves asked for thousands of deg/s
on C and B.

Each G93 move's time becomes the longer of its planned time and each axis' travel at its speed limit: moves are only
slowed down, never sped up. (The planned time already keeps the tip at Cura's speed in planar space; a stretch of
the deformation can still make the real tip faster, which nothing here limits.) G94 moves that move (the "F20000"
travel down to a re-entry point, and the first retract of a print) become G93 moves timed by their mm/min feed and
the axis limits, so every motion has a definite time.
"""
import re

import numpy as np

from .geometry import nozzle_tip

_WORD = re.compile(r"([A-Z])(-?\d+(?:\.\d*)?|-?\.\d+)")


def _fmt_f(f):
    return f"{f:.4f}" if f >= 1 else f"{f:.6g}"


def apply(path, max_c, max_b, max_x, max_z, nozzle_offset=42.0):
    """Rewrite the feeds in path in place. Speeds in deg/s (C, B) and mm/s (X, Z). Returns the share of moves slowed.
    The print-time estimate is s4/print_time.py (it also counts acceleration, like R-Theta Sim)."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    lim = {"C": max_c * 60.0, "B": max_b * 60.0, "X": max_x * 60.0, "Z": max_z * 60.0}  # per minute
    pos = {"C": 0.0, "X": 0.0, "Z": 20.0, "B": 0.0}
    out = []
    mode = 94
    started = False
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
            # a G94 / move F<mm/min> / G93 triple: the writer's feedless move (F20000), or a retract emitted at a new
            # position (the first retract of a print swings C to the first point). If it moves, make it a G93 move
            # timed by its feed and the axis limits; E-only triples (retract in place) stay as they are.
            if i + 2 < n and lines[i + 2].strip() == "G93" and started and lines[i + 1].startswith(("G00 C", "G01 C")):
                w = dict(_WORD.findall(lines[i + 1].split(";", 1)[0]))
                new = {k: float(w[k]) for k in pos}
                if "F" in w and any(new[k] != pos[k] for k in pos):
                    t = _move_time(pos, new, None, float(w["F"]), lim, nozzle_offset)
                    out.append(re.sub(r" F-?[\d.eE+-]+", " F" + _fmt_f(1.0 / t), lines[i + 1].rstrip()))
                    pos = new; moves += 1
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
                moves += 1
                if t > t0 * (1 + 1e-9):
                    slowed += 1
                    line = re.sub(r" F-?[\d.eE+-]+", " F" + _fmt_f(1.0 / t), line.rstrip())
            pos = new
        out.append(line)
        i += 1
    with open(path, "w") as fh:
        fh.write("\n".join(out))
    return {"moves_slowed_pct": round(100.0 * slowed / max(moves, 1), 1)}


def _tip(p, L):
    return np.array(nozzle_tip(p["C"], p["X"], p["Z"], p["B"], L))


def _move_time(a, b, t0, feed, lim, L):
    """Minutes for the move a -> b."""
    t = t0 or 0.0
    if feed:
        t = max(t, float(np.linalg.norm(_tip(b, L) - _tip(a, L))) / feed)
    for k, v in lim.items():
        if v > 0:
            t = max(t, abs(b[k] - a[k]) / v)
    return max(t, 1e-9)
