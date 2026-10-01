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
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from s4 import quality  # noqa: E402

MODELS = ["pi 3mm", "benchy upsidedown tilted", "Squirtle", "dino", "z mount 3mm", "B axis mount 5mm"]
# regions where roughness matters most, in the part's real frame (x, y, z boxes)
REGIONS = {"pi 3mm": {"bridge_underside": ((-5, 5), (-5.5, 5.5), (9.5, 12.5))}}

# metric: (worse when, tolerance before it counts as a regression)
CHECKS = {"ungrounded_mm": ("up", 1.0), "island_mm": ("up", 0.5), "cantilever_mm": ("up", 2.0),
          "poles": ("up", 0), "along_axis": ("up", 0), "zigzag": ("up_rel", 0.10), "tilt_deg": ("down", 2.0)}


def run_model(model, out, extra, python):
    gc = os.path.join(out, f"{model}.gcode")
    args = [python, os.path.join(HERE, "s4_slice.py"), os.path.join(HERE, "input_models", f"{model}.stl"),
            "-o", gc, "--work-dir", os.path.join(out, model)] + extra
    t = time.time()
    p = subprocess.run(args, cwd=HERE, capture_output=True, text=True, errors="replace")
    log = p.stdout + p.stderr
    open(os.path.join(out, f"{model}.log"), "w", encoding="utf-8").write(log)
    if p.returncode or not os.path.exists(gc):
        return model, {"error": (log.strip().splitlines() or ["?"])[-1][:300]}
    return model, measure(model, gc, log, time.time() - t)


def measure(model, gc, log="", seconds=None):
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
    m = re.search(r"extruding along the nozzle axis[^:]*: (\d+)", log)
    res = {k: num(k) for k in ("ungrounded_mm", "island_mm", "cantilever_mm", "bridge_mm", "poles",
                               "estimated_print_time_min")}
    res.update(tilt_deg=num("deform.tilt_deg"), target_tilt_deg=num("deform.target_tilt_deg"),
               along_axis=int(m.group(1)) if m else 0, zigzag=r["count"], segments=r["segments"],
               seconds=round(seconds, 1) if seconds else None)
    for name, box in REGIONS.get(model, {}).items():
        mid = r["mid"]; inside = (r["rough"] > 30)
        for ax, (lo, hi) in enumerate(box):
            inside &= (mid[:, ax] >= lo) & (mid[:, ax] <= hi)
        res[f"zigzag_{name}"] = int(inside.sum())
    return res


def compare(res, base):
    """Regressions of res against base (same model)."""
    bad = []
    for k, (how, tol) in CHECKS.items():
        a, b = res.get(k), base.get(k)
        if a is None or b is None:
            continue
        if (how == "up" and a > b + tol) or (how == "down" and a < b - tol) or \
                (how == "up_rel" and a > b * (1 + tol) + 5):
            bad.append(f"{k} {b:g} -> {a:g}")
    for k in res:
        if k.startswith("zigzag_") and k in base and res[k] > base[k] * 1.1 + 3:
            bad.append(f"{k} {base[k]:g} -> {res[k]:g}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", action="append", default=[], help="passed to s4_slice.py --set (repeatable)")
    ap.add_argument("--cura-set", action="append", default=[])
    ap.add_argument("--models", nargs="*", default=MODELS)
    ap.add_argument("--baseline", help="a previous --out directory to compare with")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--measure-only", action="store_true", help="don't slice: re-measure the G-code already in --out")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    extra = sum((["--set", s] for s in a.set), []) + sum((["--cura-set", s] for s in a.cura_set), [])
    results = {}
    if a.measure_only:
        for mdl in a.models:
            gc = os.path.join(a.out, f"{mdl}.gcode")
            log = open(os.path.join(a.out, f"{mdl}.log"), errors="replace").read() if os.path.exists(gc[:-6] + ".log") else ""
            results[mdl] = measure(mdl, gc, log)
    else:
        with cf.ThreadPoolExecutor(a.jobs) as ex:
            for mdl, res in ex.map(lambda m: run_model(m, a.out, extra, sys.executable), a.models):
                results[mdl] = res
    json.dump({"settings": a.set, "cura_settings": a.cura_set, "results": results},
              open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    base = json.load(open(os.path.join(a.baseline, "summary.json")))["results"] if a.baseline else {}
    cols = ["ungrounded_mm", "island_mm", "cantilever_mm", "bridge_mm", "poles", "along_axis", "tilt_deg", "zigzag",
            "estimated_print_time_min"]
    print(f"{'model':26s}" + "".join(f"{c.replace('_mm', '').replace('estimated_print_time_min', 'minutes'):>12s}" for c in cols)
          + "  region zigzag")
    regressions = 0
    for mdl in a.models:
        r = results[mdl]
        if "error" in r:
            print(f"{mdl:26s} ERROR {r['error']}"); regressions += 1; continue
        reg = " ".join(f"{k[7:]}={v}" for k, v in r.items() if k.startswith("zigzag_"))
        print(f"{mdl:26s}" + "".join(f"{('-' if r.get(c) is None else f'{r[c]:g}'):>12s}" for c in cols) + f"  {reg}")
        if mdl in base and "error" not in base[mdl]:
            bad = compare(r, base[mdl])
            regressions += bool(bad)
            print(f"{'':26s}  {'REGRESSION: ' + '; '.join(bad) if bad else 'ok vs baseline'}")
    print(f"settings: {a.set} {a.cura_set}  ->  {a.out}")
    return 1 if regressions else 0


if __name__ == "__main__":
    sys.exit(main())
