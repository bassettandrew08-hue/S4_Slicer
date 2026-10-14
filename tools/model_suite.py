"""
Slice every test model and compare quality against a baseline run, so a fix for one model can't quietly hurt another.

    venv\\Scripts\\python tools\\model_suite.py --out build/suite/base
    venv\\Scripts\\python tools\\model_suite.py --out build/suite/try1 --set SEG_SIZE=0.3 --baseline build/suite/base

Per model it reports structure (ungrounded mm, floating island / cantilever / bridge mm, poles, extrusion along the
nozzle axis), tilt, surface roughness (R-Theta Sim's path-roughness count, s4/quality.py) and the print-time estimate.
With --baseline it flags regressions. Settings come from params/<model>.json as usual, plus any --set / --cura-set.
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402

from s4 import quality  # noqa: E402
from tools import layer_lines, shape_metrics  # noqa: E402

MODELS = ["pi 3mm", "benchy upsidedown tilted", "Squirtle", "dino", "z mount 3mm", "B axis mount 5mm"]
# regions where roughness matters most, in the part's real frame (x, y, z boxes)
REGIONS = {"pi 3mm": {"bridge_underside": ((-5, 5), (-5.5, 5.5), (9.5, 12.5))}}

# the per-model metrics, in table-column order. Most are read from the run's "; s4:" header lines; tilt_deg,
# along_axis and zigzag are computed in measure()
METRICS = ("ungrounded_mm", "ungrounded_wall_mm", "island_wall_mm", "cantilever_wall_mm", "island_mm",
           "cantilever_mm", "bridge_mm", "poles", "along_axis", "tilt_deg", "zigzag", "estimated_print_time_min",
           "outside_mm", "surface_rough_pct", "asym_mm", "waviness")  # outside..asym: tools/shape_metrics.py

# metric: (worse when, tolerance before it counts as a regression)
# structure is judged on walls/skin; sparse infill floats a little even in a flat print (it is still reported)
CHECKS = {"ungrounded_mm": ("up_rel", 0.5), "ungrounded_wall_mm": ("up", 1.0), "island_wall_mm": ("up", 0.5),
          "cantilever_wall_mm": ("up", 2.0), "poles": ("up", 0), "along_axis": ("up", 0), "zigzag": ("up_rel", 0.10), "waviness": ("up_rel", 0.10),
          "tilt_deg": ("down", 2.0)}


def run_model(model, out, extra, python, sliced_from=None):
    gc = os.path.join(out, f"{model}.gcode")
    args = [python, os.path.join(HERE, "s4_slice.py"), os.path.join(HERE, "input_models", f"{model}.stl"),
            "-o", gc, "--work-dir", os.path.join(out, model)] + extra
    if sliced_from:  # map-only changes: reuse that run's planar G-code (CuraEngine isn't deterministic)
        os.makedirs(os.path.join(out, model), exist_ok=True)
        planar = os.path.join(out, model, f"{model}_deformed_tet.gcode")  # a copy, so measure() finds it here
        shutil.copyfile(os.path.join(sliced_from, model, f"{model}_deformed_tet.gcode"), planar)
        args += ["--sliced-gcode", planar]
    t = time.time()
    p = subprocess.run(args, cwd=HERE, capture_output=True, text=True, errors="replace")
    log = p.stdout + p.stderr
    open(os.path.join(out, f"{model}.log"), "w", encoding="utf-8").write(log)
    if p.returncode or not os.path.exists(gc):
        return model, {"error": (log.strip().splitlines() or ["?"])[-1][:300]}
    return model, measure(model, gc, log, time.time() - t)


def measure(model, gc, log=None, seconds=None):
    """Metrics of one finished run. log=None: the run's log is missing, so the [quality] counts are unknown."""
    if not os.path.exists(gc):
        return {"error": f"no G-code at {gc}"}
    hdr = {}
    with open(gc, errors="replace") as fh:
        for line in fh:
            m = re.match(r"; s4: (\S+) = (.*)", line)
            if m:
                hdr[m.group(1)] = m.group(2)
            elif line.startswith("; ---- end of slicer settings"):
                break
    num = lambda k, d=None: float(hdr[k]) if k in hdr and re.match(r"^-?[\d.]+(e-?\d+)?$", hdr[k]) else d
    r = quality.path_roughness(gc, num("NOZZLE_OFFSET", 42.0))
    m = re.search(r"extruding along the nozzle axis[^:]*: (\d+)", log or "")
    computed = {"tilt_deg": num("deform.tilt_deg"),
                "along_axis": int(m.group(1)) if m else (0 if log is not None else None),
                "zigzag": r["count"]}
    res = {k: computed[k] if k in computed else num(k) for k in METRICS}
    p0, p1 = quality.bead_segments(gc, num("NOZZLE_OFFSET", 42.0))
    res.update(shape_metrics.measure(model, os.path.join(HERE, "input_models", f"{model}.stl"), r["mid"],
                                     np.linalg.norm(p1 - p0, axis=1), r["rough"]))
    res.update(target_tilt_deg=num("deform.target_tilt_deg"), segments=r["segments"],
               seconds=round(seconds, 1) if seconds else None)
    # waviness: joints turning > 10 deg where the planar path runs straight (tools/roughness_breakdown.py); needs the
    # run's map_src.npy (written by --impl fast) and its planar G-code
    work = gc[:-len(".gcode")]
    src = os.path.join(work, "map_src.npy")
    planar = hdr.get("planar_gcode", f"{model}_deformed_tet.gcode")
    planar = planar if os.path.isabs(planar) else os.path.join(work, planar)
    res["waviness"] = None
    if os.path.exists(src) and os.path.exists(planar):
        from tools.roughness_breakdown import waviness
        res["waviness"] = waviness(gc, planar, np.load(src), num("SEG_SIZE", 0.6), num("RETRACTION_LENGTH", 1.0),
                                   num("NOZZLE_OFFSET", 42.0))[0]
    # layer lines on the visible outer surface (tools/layer_lines.py; regions there): needs the work dir
    if os.path.exists(os.path.join(work, "deformed_points.npy")) and os.path.exists(os.path.join(work, "params_used.json")):
        res.update(layer_lines.measure(model, work))
    for name, box in REGIONS.get(model, {}).items():
        mid = r["mid"]
        inside = (r["rough"] > 30)
        for ax, (lo, hi) in enumerate(box):
            inside &= (mid[:, ax] >= lo) & (mid[:, ax] <= hi)
        res[f"zigzag_{name}"] = int(inside.sum())
    return res


def compare(res, base):
    """Regressions of res against base (same model), and notes on metrics that could not be compared."""
    bad, notes = [], []
    for k, (how, tol) in CHECKS.items():
        a, b = res.get(k), base.get(k)
        if a is None or b is None:
            if (a is None) != (b is None):
                notes.append(f"{k} not compared ({'baseline' if b is None else 'this run'} lacks it)")
            continue
        if (how == "up" and a > b + tol) or (how == "down" and a < b - tol) or \
                (how == "up_rel" and a > b * (1 + tol) + 5):
            bad.append(f"{k} {b:g} -> {a:g}")
    for k in res:
        if k.startswith(("zigzag_", "layer_turns_")) and k in base and res[k] > base[k] * 1.1 + 3:
            bad.append(f"{k} {base[k]:g} -> {res[k]:g}")
    return bad, notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", action="append", default=[], help="passed to s4_slice.py --set (repeatable)")
    ap.add_argument("--cura-set", action="append", default=[])
    ap.add_argument("--models", nargs="*", default=MODELS)
    ap.add_argument("--baseline", help="a previous --out directory to compare with")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--sliced-from", help="reuse the planar G-code of this earlier --out (for map-only settings)")
    ap.add_argument("--measure-only", action="store_true", help="don't slice: re-measure the G-code already in --out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    extra = sum((["--set", s] for s in a.set), []) + sum((["--cura-set", s] for s in a.cura_set), [])
    results = {}
    if a.measure_only:
        for mdl in a.models:
            gc = os.path.join(a.out, f"{mdl}.gcode")
            log = open(gc[:-6] + ".log", errors="replace").read() if os.path.exists(gc[:-6] + ".log") else None
            results[mdl] = measure(mdl, gc, log)
    else:
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            for mdl, res in ex.map(lambda m: run_model(m, a.out, extra, sys.executable, a.sliced_from), a.models):
                results[mdl] = res
    json.dump({"settings": a.set, "cura_settings": a.cura_set, "results": results},
              open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    base = json.load(open(os.path.join(a.baseline, "summary.json")))["results"] if a.baseline else {}
    cols = METRICS
    heads = (c.replace("_mm", "").replace("estimated_print_time_min", "minutes") for c in cols)
    print(f"{'model':26s}" + "".join(f"{h:>12s}" for h in heads) + "  region zigzag")
    regressions = 0
    for mdl in a.models:
        r = results[mdl]
        if "error" in r:
            print(f"{mdl:26s} ERROR {r['error']}")
            regressions += 1
            continue
        reg = " ".join(f"{k[7:]}={v}" for k, v in r.items() if k.startswith("zigzag_"))
        reg += "".join(f" {k}={v}" for k, v in r.items() if k.startswith("layer_turns_"))
        print(f"{mdl:26s}" + "".join(f"{('-' if r.get(c) is None else f'{r[c]:g}'):>12s}" for c in cols) + f"  {reg}")
        if mdl in base and "error" not in base[mdl]:
            bad, notes = compare(r, base[mdl])
            regressions += bool(bad)
            print(f"{'':26s}  {'REGRESSION: ' + '; '.join(bad) if bad else 'ok vs baseline'}"
                  + (f"  (note: {'; '.join(notes)})" if notes else ""))
    print(f"settings: {a.set} {a.cura_set}  ->  {a.out}")
    return 1 if regressions else 0


if __name__ == "__main__":
    sys.exit(main())
