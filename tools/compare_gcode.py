"""
Compare two 4-axis S4 G-code files (C/X/Z/B/E/F) line by line within a tolerance.

Exit code 0 = PASS, 1 = FAIL. Usage:
    python tools/compare_gcode.py baseline.gcode candidate.gcode [--tol-pos 1e-3] ...

Tolerances are absolute, in the file's units (C,B: degrees; X,Z: mm; E: mm filament).
F (inverse-time feed) is compared relatively. Values are printed with 4-5 decimals,
so any tolerance below the print resolution is raised to one unit of the last digit.
"""
import argparse
import re
import sys

import numpy as np

WORD = re.compile(r"([A-Z])(-?\d+(?:\.\d*)?|-?\.\d+)")
LETTERS = "CXZBEF"


def parse(path):
    cmds, vals, raw = [], [], []
    with open(path) as fh:
        for line in fh:
            code = line.split(";", 1)[0].strip()
            if not code:
                continue
            parts = code.split(None, 1)
            cmd = parts[0]
            row = [np.nan] * len(LETTERS)
            if len(parts) > 1:
                for letter, num in WORD.findall(parts[1]):
                    if letter in LETTERS:
                        row[LETTERS.index(letter)] = float(num)
            cmds.append(cmd)
            vals.append(row)
            raw.append(code)
    return np.array(cmds), np.array(vals, dtype=float), raw


def compare(a_path, b_path, tol_pos=1e-3, tol_ang=1e-3, tol_e=1e-4, tol_f_rel=1e-4, quiet=False):
    out = []
    ca, va, ra = parse(a_path)
    cb, vb, rb = parse(b_path)
    ok = True
    out.append(f"baseline : {a_path}  ({len(ca)} lines)")
    out.append(f"candidate: {b_path}  ({len(cb)} lines)")
    n = min(len(ca), len(cb))
    if len(ca) != len(cb):
        ok = False
        out.append(f"FAIL: line count differs ({len(ca)} vs {len(cb)})")
    cmd_mismatch = np.nonzero(ca[:n] != cb[:n])[0]
    if len(cmd_mismatch):
        ok = False
        i = cmd_mismatch[0]
        out.append(f"FAIL: {len(cmd_mismatch)} command mismatches, first at code line {i}: "
                   f"{ra[i]!r} vs {rb[i]!r}")
        n = min(n, i)  # numeric comparison only up to first structural divergence
    A, B = va[:n], vb[:n]
    nan_mismatch = np.isnan(A) != np.isnan(B)
    if nan_mismatch.any():
        ok = False
        i = np.nonzero(nan_mismatch.any(axis=1))[0][0]
        out.append(f"FAIL: {nan_mismatch.any(axis=1).sum()} lines differ in which words are present, "
                   f"first at code line {i}: {ra[i]!r} vs {rb[i]!r}")
    # print resolution: C X Z B at 5 decimals, E F at 4
    res = np.array([1e-5, 1e-5, 1e-5, 1e-5, 1e-4, 1e-4])
    tol = np.array([tol_ang, tol_pos, tol_pos, tol_ang, tol_e, np.nan])
    tol = np.maximum(tol, res * 1.01)
    out.append(f"{'word':<5}{'max |diff|':>14}{'tol':>12}{'#over tol':>11}  worst line")
    with np.errstate(invalid="ignore"):
        diff = np.abs(A - B)
        for j, letter in enumerate(LETTERS):
            d = diff[:, j]
            if letter == "F":
                scale = np.maximum(np.abs(A[:, j]), 1e-12)
                d = d / scale
                t = tol_f_rel
                over = (d > t) & (np.abs(A[:, j] - B[:, j]) > res[j] * 1.01)
                t_show = f"{t:.0e} rel"
            else:
                t = tol[j]
                over = d > t
                t_show = f"{t:.1e}"
            dmax = np.nanmax(d) if np.isfinite(d).any() else 0.0
            worst = int(np.nanargmax(d)) if np.isfinite(d).any() else -1
            n_over = int(np.sum(over))
            if n_over:
                ok = False
            out.append(f"{letter:<5}{dmax:>14.3e}{t_show:>12}{n_over:>11}  {worst}")
    # global sanity numbers
    ea = np.nansum(np.where(va[:, 4] > 0, va[:, 4], 0))
    eb = np.nansum(np.where(vb[:, 4] > 0, vb[:, 4], 0))
    out.append(f"total positive E: {ea:.4f} vs {eb:.4f} (rel diff {abs(ea - eb) / max(ea, 1e-12):.2e})")
    out.append("PASS" if ok else "FAIL")
    if not quiet:
        print("\n".join(out))
    return ok, "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline")
    ap.add_argument("candidate")
    ap.add_argument("--tol-pos", type=float, default=1e-3, help="mm, X/Z")
    ap.add_argument("--tol-ang", type=float, default=1e-3, help="deg, C/B")
    ap.add_argument("--tol-e", type=float, default=1e-4, help="mm filament")
    ap.add_argument("--tol-f-rel", type=float, default=1e-4, help="relative, F")
    a = ap.parse_args()
    ok, _ = compare(a.baseline, a.candidate, a.tol_pos, a.tol_ang, a.tol_e, a.tol_f_rel)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
