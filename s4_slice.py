"""
One command: STL in, 4-axis (C/X/Z/B/E, G93 inverse-time, M83) G-code out.

    venv\\Scripts\\python s4_slice.py "input_models/pi 3mm.stl"
    venv\\Scripts\\python s4_slice.py model.stl --set MAX_OVERHANG=10 --cura-set layer_height=0.1

Stages: deform (tetgen + rotation field + deformation solve) -> CuraEngine (headless, settings
from cura_config.3mf) -> map planar G-code back to the 4-axis machine.

Per-build settings live in a profile: params/<model name>.json is used automatically if it exists
(or pass --params FILE). Create one with --init-params, check the result with --show-params.
The settings each run actually used are saved to build/<model>/params_used.json.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from s4 import pipeline  # noqa: E402
from s4 import profile as profiles  # noqa: E402
from s4.timing import TIMER  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", help="input STL")
    ap.add_argument("-o", "--out", help="output G-code (default: output_gcode/<model>.gcode)")

    g = ap.add_argument_group("settings (later ones win: defaults < profile < --set < --cura-*)")
    g.add_argument("--params", metavar="FILE", help="build profile JSON (default: params/<model>.json if it exists)")
    g.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", dest="sets",
                   help="override a deform or map setting, e.g. --set MAX_OVERHANG=10 --set PART_OFFSET=[0,10,0]")
    g.add_argument("--cura-set", action="append", default=[], metavar="KEY=VALUE", help="override a Cura setting")
    g.add_argument("--cura-config", help="Cura project (.3mf) with the base profile")
    g.add_argument("--notebook-exact", action="store_true",
                   help="reproduce the notebook exactly: its deformation solve, and none of the retraction, "
                        "extrusion-multiplier or start-prime fixes")
    g.add_argument("--show-params", action="store_true", help="print the resolved settings and exit")
    g.add_argument("--init-params", nargs="?", const="", metavar="FILE",
                   help="write the resolved settings to a profile (default: params/<model>.json) and exit")
    g.add_argument("--force", action="store_true", help="let --init-params overwrite an existing file")

    o = ap.add_argument_group("run options")
    o.add_argument("--impl", choices=["fast", "reference"], default="fast",
                   help="fast (default) or the verified line-for-line notebook port")
    o.add_argument("--work-dir", help="intermediate files (default: build/<model>/)")
    o.add_argument("--cura-engine", help="path to CuraEngine.exe (default: newest UltiMaker Cura install)")
    o.add_argument("--sliced-gcode", help="skip CuraEngine and map this planar G-code instead")
    o.add_argument("--save-gif", action="store_true", help="write the notebook's progress GIFs (reference impl only)")
    o.add_argument("--save-pickle", action="store_true", help="also pickle the deformed tet mesh like the notebook")
    o.add_argument("--no-support-check", action="store_true", help="skip the mid-air extrusion report")
    a = ap.parse_args()

    name = os.path.splitext(os.path.basename(a.model))[0]
    cura_sets = {}
    for kv in a.cura_set:
        k, sep, v = kv.partition("=")
        if not sep:
            sys.exit(f"--cura-set expects KEY=VALUE, got {kv!r}")
        cura_sets[k.strip()] = v.strip()
    try:
        prof, source = profiles.resolve(a.model, a.params, a.sets, cura_sets, a.cura_config, a.notebook_exact)
    except (profiles.ProfileError, OSError) as e:
        sys.exit(f"error: {e}")
    if prof["cura"]["set"] and not a.sliced_gcode:
        from s4.cura import check_setting_names
        try:
            check_setting_names(profiles.cura_config_path(prof), prof["cura"]["set"], a.cura_engine)
        except (ValueError, OSError) as e:
            sys.exit(f"error: {e}")

    if a.init_params is not None:
        path = a.init_params or profiles.default_profile_path(a.model)
        if os.path.exists(path) and not a.force:
            sys.exit(f"error: {path} already exists (use --force to overwrite)")
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(profiles.to_json(prof) + "\n")
        print(f"wrote {path}  (settings from: {source}). Edit it, then run without --init-params.")
        return
    if a.show_params:
        print(f"# settings from: {source}" + (f" + {len(a.sets)} --set" if a.sets else ""))
        print(profiles.to_json(prof))
        return

    print(f"[params] {source}" + (f" + --set {' '.join(a.sets)}" if a.sets else "")
          + (f" + --cura-set {' '.join(a.cura_set)}" if a.cura_set else ""))
    out = a.out or os.path.join(HERE, "output_gcode", f"{name}.gcode")
    pipeline.run(a.model, out, prof, impl=a.impl, work_dir=a.work_dir or os.path.join(HERE, "build", name),
                 cura_engine=a.cura_engine, sliced_gcode=a.sliced_gcode, save_gif=a.save_gif,
                 save_pickle=a.save_pickle, support_check=not a.no_support_check)
    print(TIMER.report(f"Stage timings: {name} ({a.impl})"))


if __name__ == "__main__":
    main()
