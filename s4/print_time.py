"""
Print-time estimate that matches R-Theta Sim (sim/r-theta-simulator.html).

A line-by-line port of the sim's parseGcode() and plan(): per-axis speed and acceleration limits, trapezoidal
moves, and junction speeds from a backward and a forward lookahead pass. Keep the two in step: a change to the
sim's planner needs the same change here (tools/check_print_time.py compares them).

The sim's quirks are copied on purpose, so both give the same number:
- a move's length is its X/Z distance in mm; if only B/C move it is their distance in degrees; if only E moves, the
  E distance. Speeds, accelerations and the corner-speed floor apply in those units.
- junctions only carry speed between two moves of the same kind (X/Z, B/C only, E only); anything else stops.
- F is sticky from 1500 mm/min; the machine starts at X Z B C E = 0 and G28 moves to the home position.
"""
import math
import re

import numpy as np

try:
    import numba
    _njit = numba.njit(cache=True)
except ImportError:  # pragma: no cover
    def _njit(f):
        return f

# the machine values the planner uses (pipeline map setting names; the sim's vX.., aX.., jerk, home*).
# Same defaults as s4/profile.py MAP_DEFAULTS and the sim's DEFAULTS
MACHINE_DEFAULTS = dict(
    MAX_SPEED_X=150, MAX_SPEED_Z=50, MAX_SPEED_B=180, MAX_SPEED_C=360, MAX_SPEED_E=60,
    MAX_ACCEL_X=2000, MAX_ACCEL_Z=200, MAX_ACCEL_B=5000, MAX_ACCEL_C=3000, MAX_ACCEL_E=3000,
    CORNER_SPEED=8, HOME_X=140, HOME_Z=200, HOME_B=0, HOME_SPEED=40,
)

_RE = re.compile(r"([A-Z])\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)", re.ASCII)  # JS \d is ASCII only
_PAREN = re.compile(r"\(.*?\)")
_AX = {"X": 0, "Z": 1, "B": 2, "C": 3}


def parse(lines, m):
    """The sim's parseGcode (feedMode 'linear', g0Inverse 'inverse'). Returns arrays like the sim's P."""
    st = [0.0, 0.0, 0.0, 0.0, 0.0]
    pts = [st[:]]
    typ, feed, dwell, inv = [], [], [], []
    abs_xyz, abs_e, F, last_motion, inv_mode = True, True, 1500.0, 1, False
    home = (float(m["HOME_X"]), float(m["HOME_Z"]), float(m["HOME_B"]))

    def push(t, f, dw, s, iv):
        typ.append(t); feed.append(f); dwell.append(dw); inv.append(iv)
        st[:] = s; pts.append(s[:])

    for raw in lines:
        s = raw
        sc = s.find(";")
        if sc >= 0:
            s = s[:sc]
        if "(" in s:
            s = _PAREN.sub("", s)
        s = s.strip()
        if not s:
            continue
        s = s.upper()
        g = mc = -1.0
        w = {}
        has_axis = False
        for L, v in _RE.findall(s):
            v = float(v)
            if L == "G" and g < 0:
                g = v
            elif L == "M" and mc < 0:
                mc = v
            else:
                w[L] = v
                if L in _AX or L == "E" or L == "Y":
                    has_axis = True
        if g == 0 or g == 1 or (g < 0 and mc < 0 and has_axis):
            gg = g if g >= 0 else last_motion
            last_motion = gg
            if "F" in w:
                F = w["F"]
            t = st[:]
            for k, j in _AX.items():
                if k in w:
                    t[j] = w[k] if abs_xyz else st[j] + w[k]
            if "E" in w:
                t[4] = w["E"] if abs_e else st[4] + w["E"]
            if t != st:
                push(0 if gg == 0 else 1, F, 0.0, t, 1 if inv_mode else 0)
        elif g == 28:
            t = st[:]
            any_ax = any(k in w for k in ("X", "Z", "B", "C"))
            if not any_ax or "X" in w:
                t[0] = home[0]
            if not any_ax or "Z" in w:
                t[1] = home[1]
            if not any_ax or "B" in w:
                t[2] = home[2]
            if "C" in w:
                t[3] = 0.0
            push(2, F, 0.0, t, 0)
        elif g == 92:
            t = st[:]
            for k, j in _AX.items():
                if k in w:
                    t[j] = w[k]
            if "E" in w:
                t[4] = w["E"]
            push(3, F, 0.0, t, 0)
        elif g == 4:
            secs = w["S"] if "S" in w else (w["P"] / 1000 if "P" in w else 0.0)
            if secs > 0:
                push(4, F, secs, st[:], 0)
        elif g == 93:
            inv_mode = True
        elif g == 94:
            inv_mode = False
        elif g == 90:
            abs_xyz = True
        elif g == 91:
            abs_xyz = False
        elif mc == 82:
            abs_e = True
        elif mc == 83:
            abs_e = False
    # the sim stores feed and dwell as Float32Array
    return (np.array(typ, np.uint8), np.array(feed, np.float32).astype(np.float64),
            np.array(dwell, np.float32).astype(np.float64), np.array(inv, np.uint8), np.array(pts, np.float64))


@_njit
def _plan(typ, feed, dwell, inv, c, vmax, amax, vE, aE, jerk, home_speed):
    n = len(typ)
    L = np.zeros(n); vc = np.zeros(n); acc = np.zeros(n); kind = np.zeros(n, np.int8); U = np.zeros((n, 4))
    d = np.zeros(4)
    for i in range(n):
        t = typ[i]
        if t == 3 or t == 4:
            continue
        for k in range(4):
            d[k] = c[i + 1, k] - c[i, k]
        dE = c[i + 1, 4] - c[i, 4]
        lin = math.hypot(d[0], d[1]); rot = math.hypot(d[2], d[3])
        if lin > 1e-9:
            ln = lin; kd = 1
        elif rot > 1e-9:
            ln = rot; kd = 2
        elif abs(dE) > 1e-9:
            ln = abs(dE); kd = 3
        else:
            continue
        if t == 2:
            v = home_speed
        elif inv[i]:
            v = max(ln * feed[i] / 60.0, 1e-4)  # G93: finish the move in 1/F minutes
        else:
            v = max(feed[i] / 60.0, 0.01)
        T = ln / v; a = np.inf
        for k in range(4):
            ad = abs(d[k])
            if ad > 1e-12:
                T = max(T, ad / vmax[k]); a = min(a, amax[k] * ln / ad)
        if abs(dE) > 1e-12:
            T = max(T, abs(dE) / vE); a = min(a, aE * ln / abs(dE))
        L[i] = ln; vc[i] = ln / T; acc[i] = a; kind[i] = kd
        nrm = 0.0
        for k in range(4):
            U[i, k] = d[k] / vmax[k]; nrm += U[i, k] ** 2
        nrm = math.sqrt(nrm)
        if nrm == 0.0:
            nrm = 1.0
        for k in range(4):
            U[i, k] /= nrm
    J = np.zeros(n + 1)
    for i in range(1, n):
        if kind[i] > 0 and kind[i] == kind[i - 1]:
            cs = 0.0
            for k in range(4):
                cs += U[i, k] * U[i - 1, k]
            vmin = min(vc[i], vc[i - 1])
            J[i] = max(vmin * max(0.0, cs), min(jerk, vmin))
    for i in range(n - 1, -1, -1):
        if L[i] > 0:
            J[i] = min(J[i], math.sqrt(J[i + 1] ** 2 + 2 * acc[i] * L[i]))
    for i in range(n):
        if L[i] > 0:
            J[i + 1] = min(J[i + 1], math.sqrt(J[i] ** 2 + 2 * acc[i] * L[i]))
    time = 0.0
    for i in range(n):
        if typ[i] == 4:
            time += dwell[i]
            continue
        if L[i] <= 0:
            continue
        a = acc[i]; v0 = J[i]; v1 = J[i + 1]; ln = L[i]; vp = vc[i]
        fin = np.isfinite(a)
        dA = (vp * vp - v0 * v0) / (2 * a); dD = (vp * vp - v1 * v1) / (2 * a)
        if dA + dD > ln:
            vp = math.sqrt((2 * a * ln + v0 * v0 + v1 * v1) / 2)
            dA = (vp * vp - v0 * v0) / (2 * a); dD = (vp * vp - v1 * v1) / (2 * a)
        tA = (vp - v0) / a if fin else 0.0
        tD = (vp - v1) / a if fin else 0.0
        tC = max(0.0, (ln - dA - dD) / vp)
        time += tA + tC + tD
    return time


def estimate(path_or_lines, machine=None):
    """Seconds R-Theta Sim would take for this G-code with these machine limits (pipeline map setting names)."""
    m = {**MACHINE_DEFAULTS, **{k: v for k, v in (machine or {}).items() if k in MACHINE_DEFAULTS}}
    for k in ("MAX_SPEED_X", "MAX_SPEED_Z", "MAX_SPEED_B", "MAX_SPEED_C", "MAX_SPEED_E", "MAX_ACCEL_X", "MAX_ACCEL_Z",
              "MAX_ACCEL_B", "MAX_ACCEL_C", "MAX_ACCEL_E", "HOME_SPEED"):
        if not float(m[k]) > 0:
            raise ValueError(f"{k} must be > 0 (got {m[k]!r})")
    if isinstance(path_or_lines, str):
        with open(path_or_lines, encoding="utf-8", errors="replace", newline="") as fh:
            lines = re.split(r"\r?\n", fh.read())  # the sim splits on /\r?\n/
    else:
        lines = path_or_lines
    typ, feed, dwell, inv, c = parse(lines, m)
    if not len(typ):
        return 0.0
    vmax = np.array([m["MAX_SPEED_X"], m["MAX_SPEED_Z"], m["MAX_SPEED_B"], m["MAX_SPEED_C"]], np.float64)
    amax = np.array([m["MAX_ACCEL_X"], m["MAX_ACCEL_Z"], m["MAX_ACCEL_B"], m["MAX_ACCEL_C"]], np.float64)
    return float(_plan(typ, feed, dwell, inv, c, vmax, amax, float(m["MAX_SPEED_E"]), float(m["MAX_ACCEL_E"]),
                       float(m["CORNER_SPEED"]), float(m["HOME_SPEED"])))
