"""STL -> deform -> CuraEngine -> 4-axis map -> G-code, with per-stage timings."""
import json
import os
import pickle
import time

import numpy as np

from .timing import TIMER
from . import profile as profiles


def run(model_path, out_gcode, profile=None, impl="fast", work_dir=None, cura_engine=None, sliced_gcode=None,
        save_gif=False, save_pickle=False, support_check=True, params=None, log=print):
    """
    profile: build profile from s4.profile.resolve() (deform / map / cura settings). None = defaults.
    params: shortcut for a partial profile dict (e.g. {"deform": {...}} or a legacy flat params dict).
    impl: "fast" (default) or "reference" (the verified notebook port, slow).
    sliced_gcode: use this planar G-code instead of running CuraEngine (e.g. for A/B checks).
    support_check: report extrusion that would be printed in mid-air (s4/support_check.py).
    The settings actually used are written to <work_dir>/params_used.json (re-usable with --params).
    Returns dict with paths, stats and timings.
    """
    prof = profile or profiles.defaults()
    if params:
        prof = profiles.merge(prof, params, where="params")
    p = profiles.deform_params(prof)
    name = os.path.splitext(os.path.basename(model_path))[0]
    work_dir = work_dir or os.path.join("build", name)
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(out_gcode)), exist_ok=True)
    with open(os.path.join(work_dir, "params_used.json"), "w", encoding="utf-8") as fh:
        fh.write(profiles.to_json(prof))
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
                overrides = {k: str(v) if not isinstance(v, str) else v for k, v in prof["cura"]["set"].items()}
                if prof["cura"]["strip_start_prime"] and "machine_start_gcode" not in overrides:
                    # Cura's start code primes (G1 F200 E3) at the park position, which the mapper would place
                    # inside the part as a floating blob; the S4 header already primes at home (G1 E10).
                    overrides["machine_start_gcode"] = "G28 ; home"
                cura_info = slice_stl(stl_path, planar_path, profiles.cura_config_path(prof), overrides, cura_engine, log=log)
        else:
            log(f"[slice] using existing planar G-code {sliced_gcode}")
        retraction = 1.0
        if cura_info:
            retraction = float(cura_info["extruder"].get("retraction_amount", cura_info["global"].get("retraction_amount")))

        # ---- 3. map back to 4 axes
        with TIMER("3. map to 4-axis G-code"):
            mp = dict(prof["map"])
            if mp["RETRACTION_LENGTH"] is None:
                mp["RETRACTION_LENGTH"] = retraction
            retraction = mp["RETRACTION_LENGTH"]
            if impl == "reference":
                from . import reference as ref
                pts, stats = ref.map_gcode(input_tet, deformed, planar_path, mp)
                with TIMER("write G-code"):
                    ref.write_gcode(pts, out_gcode, NOZZLE_OFFSET=mp["NOZZLE_OFFSET"])
            else:
                from . import fast_map, meshio_s4
                d_points = np.asarray(deformed.points)
                d_centers = np.asarray(deformed.cell_centers().points)
                pts, stats = fast_map.map_gcode(in_cells, in_points, in_centers, d_points, d_centers,
                                                meshio_s4.make_find_cells(deformed), planar_path, mp)
                with TIMER("write G-code"):
                    fast_map.write_gcode(pts, out_gcode, NOZZLE_OFFSET=mp["NOZZLE_OFFSET"])
        if support_check:
            with TIMER("4. support check"):
                from . import support_check as sc
                from .params import expand_iterations
                support = sc.check(model_path, np.asarray(deformed.points), planar_path,
                                   part_offset=expand_iterations(p)[0], retraction_length=retraction,
                                   seg_size=mp["SEG_SIZE"])
            stats["floating_points"] = support["floating_points"]
            log(sc.format_report(support))
    from . import support_check as sc2
    poles = sc2.vertical_extrusion(out_gcode, nozzle_offset=mp["NOZZLE_OFFSET"])
    stats["poles"] = len(poles[0])
    log(sc2.format_vertical(poles))
    total = time.perf_counter() - t0
    log(f"[map] {stats}")
    log(f"[done] {out_gcode}  ({total:.1f} s)")
    result = {"out": out_gcode, "stl": stl_path, "planar": planar_path, "stats": stats, "total_s": total,
              "timings": TIMER.as_dict()}
    with open(os.path.join(work_dir, f"timings_{impl}.json"), "w") as fh:
        json.dump(result, fh, indent=1, default=str)
    return result
