"""
Run the faithful notebook port (s4/reference.py) and save artifacts + timings.

  deform : cells 2-11  -> <out>/deformed_points.npy, rotation_field.npy, deformed.stl, deformed_tet.pkl
  map    : cells 14-18 -> <out>/output.gcode   (needs --sliced G-code and a deformed mesh)

The deformed mesh for `map` comes from --deformed-pkl (e.g. the notebook's
pickle_files/deformed_<model>.pkl) or from <out>/deformed_tet.pkl.
"""
import argparse
import json
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402

from s4 import reference as ref  # noqa: E402
from s4.params import DEFAULT_PARAMS  # noqa: E402
from s4.timing import TIMER  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["deform", "map", "all"])
    ap.add_argument("--model", required=True, help="path to input STL")
    ap.add_argument("--out", required=True, help="artifact directory")
    ap.add_argument("--sliced", help="planar G-code of the deformed STL (for map)")
    ap.add_argument("--deformed-pkl", help="pickled deformed_tet to map with")
    ap.add_argument("--save-gif", action="store_true")
    ap.add_argument("--params", help="JSON file overriding DEFAULT_PARAMS")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    params = dict(DEFAULT_PARAMS)
    if args.params:
        params.update(json.load(open(args.params)))
    name = os.path.splitext(os.path.basename(args.model))[0]
    t_start = time.perf_counter()

    if args.stage in ("deform", "all"):
        with TIMER("DEFORM total"):
            input_tet, deformed_tet, rotation_field = ref.deform(args.model, name, params, save_gif=args.save_gif)
        with TIMER("save STL + pickle"):
            ref.save_deformed_stl(deformed_tet, os.path.join(args.out, "deformed.stl"))
            with open(os.path.join(args.out, "deformed_tet.pkl"), "wb") as f:
                pickle.dump(deformed_tet, f)
        np.save(os.path.join(args.out, "rotation_field.npy"), rotation_field)
        np.save(os.path.join(args.out, "deformed_points.npy"), np.asarray(deformed_tet.points))
        print(f"cells={input_tet.number_of_cells} points={input_tet.number_of_points}")

    if args.stage in ("map", "all"):
        if not args.sliced:
            sys.exit("--sliced is required for map")
        if args.stage == "map":
            with TIMER("cell 2: mesh setup (for input_tet)"):
                input_tet = ref.run_cell2(args.model, np.asarray(params["PART_OFFSET"], dtype=float), name,
                                          params.get("SURFACE_SIMPLIFY_ERROR", 0.0))
        pkl = args.deformed_pkl or os.path.join(args.out, "deformed_tet.pkl")
        with TIMER("load deformed pickle"):
            deformed_tet = pickle.load(open(pkl, "rb"))
        with TIMER("MAP total"):
            new_points, stats = ref.map_gcode(input_tet, deformed_tet, args.sliced)
            with TIMER("write 4-axis gcode"):
                ref.write_gcode(new_points, os.path.join(args.out, "output.gcode"))
        print(json.dumps(stats))
        json.dump(stats, open(os.path.join(args.out, "map_stats.json"), "w"), indent=1)

    total = time.perf_counter() - t_start
    print(TIMER.report(f"Reference timings ({name}, {args.stage}) total {total:.1f}s"))
    json.dump({"total": total, "stages": TIMER.as_dict()},
              open(os.path.join(args.out, f"timings_{args.stage}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
