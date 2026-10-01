"""
Regression gate: the fast pipeline must reproduce the reference (verified notebook port).

    venv\\Scripts\\python tools\\check_equivalence.py "input_models/pi 3mm.stl"

1. runs the reference pipeline (deform + CuraEngine + map)        -> build/check/<model>/reference/
2. runs the fast pipeline on the same model, mapping the *same* planar G-code
3. compares deformed vertices, the deformed STL, and the final 4-axis G-code (tools/compare_gcode.py)

Exit code 0 = PASS. Use --reference-dir to reuse an existing reference run.
"""
import argparse
import filecmp
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from tools.compare_gcode import compare  # noqa: E402

PY = sys.executable


def run(args):
    print(">", " ".join(f'"{a}"' if " " in a else a for a in args), flush=True)
    r = subprocess.run(args, cwd=HERE)
    if r.returncode:
        sys.exit(f"command failed ({r.returncode})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--reference-dir", help="reuse this reference work dir instead of re-running the slow reference")
    ap.add_argument("--impl", default="fast")
    ap.add_argument("--tol-pos", type=float, default=1e-3)
    ap.add_argument("--tol-ang", type=float, default=1e-3)
    ap.add_argument("--tol-e", type=float, default=1e-4)
    ap.add_argument("--params")
    ap.add_argument("--notebook-exact", action="store_true", help="pass --notebook-exact to both runs")
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[], help="extra args for the candidate run")
    a = ap.parse_args()

    name = os.path.splitext(os.path.basename(a.model))[0]
    base = os.path.join(HERE, "build", "check", name)
    ref_dir = a.reference_dir or os.path.join(base, "reference")
    cand_dir = os.path.join(base, a.impl)
    ref_out = os.path.join(ref_dir, f"{name}.gcode")
    cand_out = os.path.join(cand_dir, f"{name}.gcode")
    extra = ["--params", a.params] if a.params else []
    if a.notebook_exact:
        extra.append("--notebook-exact")

    if not a.reference_dir or not os.path.exists(ref_out):
        run([PY, "s4_slice.py", a.model, "--impl", "reference", "--work-dir", ref_dir, "-o", ref_out] + extra)
    planar = os.path.join(ref_dir, f"{name}_deformed_tet.gcode")
    run([PY, "s4_slice.py", a.model, "--impl", a.impl, "--work-dir", cand_dir, "-o", cand_out,
         "--sliced-gcode", planar] + extra + a.extra)

    ok = True
    pa = np.load(os.path.join(ref_dir, "deformed_points.npy"))
    pb = np.load(os.path.join(cand_dir, "deformed_points.npy"))
    dmax = float(np.abs(pa - pb).max()) if pa.shape == pb.shape else float("inf")
    print(f"deformed vertices: max |diff| = {dmax:.3e} mm  (bit-identical: {np.array_equal(pa, pb)})")
    ok &= dmax <= a.tol_pos
    same_stl = filecmp.cmp(os.path.join(ref_dir, f"{name}_deformed_tet.stl"),
                           os.path.join(cand_dir, f"{name}_deformed_tet.stl"), shallow=False)
    print(f"deformed STL byte-identical: {same_stl}"
          + ("" if same_stl else "  (planar G-code from CuraEngine may then differ; the map check below uses the reference slice)"))
    g_ok, _ = compare(ref_out, cand_out, a.tol_pos, a.tol_ang, a.tol_e)
    ok &= g_ok
    same_gcode = filecmp.cmp(ref_out, cand_out, shallow=False)
    print(f"final G-code byte-identical: {same_gcode}")
    t = {}
    for label, d in (("reference", ref_dir), (a.impl, cand_dir)):
        for fn in os.listdir(d):
            if fn.startswith("timings_") and fn.endswith(".json"):
                t[label] = json.load(open(os.path.join(d, fn)))["total_s"]
    if len(t) == 2:
        print(f"wall time: reference {t['reference']:.1f} s, {a.impl} {t[a.impl]:.1f} s")
    print("EQUIVALENCE", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
