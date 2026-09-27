"""
One command: STL in, 4-axis (C/X/Z/B/E, G93 inverse-time, M83) G-code out.

    venv\\Scripts\\python s4_slice.py "input_models/pi 3mm.stl"
    venv\\Scripts\\python s4_slice.py model.stl -o out.gcode --cura-set layer_height=0.1

Stages: deform (tetgen + rotation field + deformation solve) -> CuraEngine (headless, settings
from cura_config.3mf) -> map planar G-code back to the 4-axis machine.
Deformation parameters default to the values in main.ipynb; override with --params file.json.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from s4 import pipeline  # noqa: E402
from s4.timing import TIMER  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", help="input STL")
    ap.add_argument("-o", "--out", help="output G-code (default: output_gcode/<model>.gcode)")
    ap.add_argument("--impl", choices=["fast", "reference"], default="fast",
                    help="fast (default) or the verified line-for-line notebook port")
    ap.add_argument("--work-dir", help="intermediate files (default: build/<model>/)")
    ap.add_argument("--cura-config", default=os.path.join(HERE, "cura_config.3mf"), help="Cura project with the profile")
    ap.add_argument("--cura-engine", help="path to CuraEngine.exe (default: newest UltiMaker Cura install)")
    ap.add_argument("--cura-set", action="append", default=[], metavar="KEY=VALUE", help="override a Cura setting")
    ap.add_argument("--sliced-gcode", help="skip CuraEngine and map this planar G-code instead")
    ap.add_argument("--params", help="JSON with deformation parameter overrides (see s4/params.py)")
    ap.add_argument("--save-gif", action="store_true", help="write the notebook's progress GIFs (reference impl only)")
    ap.add_argument("--save-pickle", action="store_true", help="also pickle the deformed tet mesh like the notebook")
    ap.add_argument("--notebook-exact", action="store_true",
                    help="reproduce the notebook exactly: extrude the 1 mm (un)retraction during the travel "
                         "lift/plunge, the per-tet extrusion multiplier and Cura's start-code prime")
    ap.add_argument("--no-support-check", action="store_true", help="skip the mid-air extrusion report")
    a = ap.parse_args()

    name = os.path.splitext(os.path.basename(a.model))[0]
    out = a.out or os.path.join(HERE, "output_gcode", f"{name}.gcode")
    overrides = {}
    for kv in a.cura_set:
        k, sep, v = kv.partition("=")
        if not sep:
            sys.exit(f"--cura-set expects KEY=VALUE, got {kv!r}")
        overrides[k.strip()] = v.strip()
    params = json.load(open(a.params)) if a.params else None

    pipeline.run(a.model, out, impl=a.impl, work_dir=a.work_dir or os.path.join(HERE, "build", name),
                 cura_config=a.cura_config, cura_overrides=overrides, cura_engine=a.cura_engine, params=params,
                 sliced_gcode=a.sliced_gcode, save_gif=a.save_gif, save_pickle=a.save_pickle, notebook_exact=a.notebook_exact,
                 support_check=not a.no_support_check)
    print(TIMER.report(f"Stage timings: {name} ({a.impl})"))


if __name__ == "__main__":
    main()
