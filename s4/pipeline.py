"""STL -> deform -> CuraEngine -> 4-axis map -> G-code, with per-stage timings."""
import json
import os
import pickle
import time

import numpy as np

from .timing import TIMER
from .params import DEFAULT_PARAMS


def run(model_path, out_gcode, impl="fast", work_dir=None, cura_config="cura_config.3mf", cura_overrides=None,
        cura_engine=None, params=None, sliced_gcode=None, save_gif=False, save_pickle=False, log=print):
    """
    impl: "fast" (default) or "reference" (the verified notebook port, slow).
    sliced_gcode: use this planar G-code instead of running CuraEngine (e.g. for A/B checks).
    Returns dict with paths, stats and timings.
    """
    p = dict(DEFAULT_PARAMS)
    p.update(params or {})
    name = os.path.splitext(os.path.basename(model_path))[0]
    work_dir = work_dir or os.path.join("build", name)
    os.makedirs(work_dir, exist_ok=True)
    stl_path = os.path.join(work_dir, f"{name}_deformed_tet.stl")
    planar_path = sliced_gcode or os.path.join(work_dir, f"{name}_deformed_tet.gcode")
    t0 = time.perf_counter()

    with TIMER(f"TOTAL ({impl})"):
        # ---- 1. deform
        with TIMER("1. deform"):
            if impl == "reference":
                from . import reference as ref
                input_tet, deformed, rf = ref.deform(model_path, name, p, save_gif=save_gif, verbose=0)
                in_cells = input_tet.field_data["cells"]
                in_points = np.asarray(input_tet.points)
                in_centers = np.asarray(input_tet.cell_data["cell_center"])
            else:
                if save_gif:
                    log("[deform] --save-gif is only supported with --impl reference; ignoring")
                from . import fast_deform
                ctx, deformed, rf = fast_deform.deform(model_path, p)
                input_tet = ctx.tet
                in_cells = ctx.cells
                in_points = np.asarray(input_tet.points)
                in_centers = np.asarray(input_tet.cell_data["cell_center"])
        with TIMER("   save deformed STL"):
            deformed.extract_surface(algorithm='dataset_surface').save(stl_path)
            np.save(os.path.join(work_dir, "deformed_points.npy"), np.asarray(deformed.points))
            if save_pickle:
                with open(os.path.join(work_dir, f"deformed_{name}.pkl"), "wb") as f:
                    pickle.dump(deformed, f)
        log(f"[deform] {input_tet.number_of_cells} tets, {input_tet.number_of_points} vertices -> {stl_path}")

        # ---- 2. slice
        cura_info = None
        if sliced_gcode is None:
            with TIMER("2. slice (CuraEngine)"):
                from .cura import slice_stl
                cura_info = slice_stl(stl_path, planar_path, cura_config, cura_overrides, cura_engine, log=log)
        else:
            log(f"[slice] using existing planar G-code {sliced_gcode}")
        retraction = 1.0
        if cura_info:
            retraction = float(cura_info["extruder"].get("retraction_amount", cura_info["global"].get("retraction_amount")))

        # ---- 3. map back to 4 axes
        with TIMER("3. map to 4-axis G-code"):
            mp = {"RETRACTION_LENGTH": retraction}
            if impl == "reference":
                from . import reference as ref
                pts, stats = ref.map_gcode(input_tet, deformed, planar_path, mp)
                with TIMER("write G-code"):
                    ref.write_gcode(pts, out_gcode)
            else:
                from . import fast_map, meshio_s4
                d_points = np.asarray(deformed.points)
                d_centers = np.asarray(deformed.cell_centers().points)
                pts, stats = fast_map.map_gcode(in_cells, in_points, in_centers, d_points, d_centers,
                                                meshio_s4.make_find_cells(deformed), planar_path, mp)
                with TIMER("write G-code"):
                    fast_map.write_gcode(pts, out_gcode)
    total = time.perf_counter() - t0
    log(f"[map] {stats}")
    log(f"[done] {out_gcode}  ({total:.1f} s)")
    result = {"out": out_gcode, "stl": stl_path, "planar": planar_path, "stats": stats, "total_s": total,
              "timings": TIMER.as_dict()}
    with open(os.path.join(work_dir, f"timings_{impl}.json"), "w") as fh:
        json.dump(result, fh, indent=1, default=str)
    return result
