"""STL -> deform -> CuraEngine -> 4-axis map -> G-code, with per-stage timings."""
import json
import os
import pickle
import time

import numpy as np

from .timing import TIMER
from . import profile as profiles


def planar_retraction(planar_path, default=1.0):
    """Retraction length of a planar G-code: the most common E of its E-only retract lines (Cura: `G1 F.. E-1`)."""
    import collections
    import re
    pat = re.compile(r"^G1\s+(?:F[\d.]+\s+)?E-([\d.]+)\s*(?:;.*)?$")
    counts = collections.Counter()
    with open(planar_path, errors="replace") as fh:
        for line in fh:
            m = pat.match(line)
            if m:
                counts[float(m.group(1))] += 1
    return counts.most_common(1)[0][0] if counts else default


def run(model_path, out_gcode, profile=None, impl="fast", work_dir=None, cura_engine=None, sliced_gcode=None,
        save_gif=False, save_pickle=False, support_check=True, params=None, log=print):
    """
    STL -> 4-axis G-code at out_gcode. Stages:
      1. deform: tetgen mesh, rotation field, deformed mesh (fast_deform / island_free, or reference.deform)
         -> <work_dir>/<model>_deformed_tet.stl and deformed_points.npy
      2. slice: CuraEngine on the deformed STL (skipped with sliced_gcode)
      3. map: planar G-code back to the 4-axis machine (fast_map, or reference.map_gcode), then the axis speed
         limits (feed_limits, if LIMIT_AXIS_SPEEDS)
      4. support check (support_check), then the poles / along-the-nozzle-axis check (quality) and the
         print-time estimate (print_time)
      5. settings header for R-Theta Sim (sim_header) prepended to the G-code; timings_<impl>.json
    The fast and reference branches must give identical results (tools/check_equivalence.py).

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
    overrides = {k: str(v) if not isinstance(v, str) else v for k, v in prof["cura"]["set"].items()}
    if prof["cura"]["strip_start_prime"] and "machine_start_gcode" not in overrides:
        # Cura's start code primes (G1 F200 E3) at the park position, which the mapper would place
        # inside the part as a floating blob; the S4 header already primes at home (G1 E10).
        overrides["machine_start_gcode"] = "G28 ; home"
    name = os.path.splitext(os.path.basename(model_path))[0]
    work_dir = work_dir or os.path.join(profiles.HERE, "build", name)  # profiles.HERE: the repo root
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(out_gcode)), exist_ok=True)
    with open(os.path.join(work_dir, "params_used.json"), "w", encoding="utf-8") as fh:
        fh.write(profiles.to_json(prof))
    stl_path = os.path.join(work_dir, f"{name}_deformed_tet.stl")
    planar_path = sliced_gcode or os.path.join(work_dir, f"{name}_deformed_tet.gcode")
    t0 = time.perf_counter()
    from . import support_check as sc  # heavy imports (open3d, pyvista): only when a run starts
    from . import quality

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
                cura_info = slice_stl(stl_path, planar_path, profiles.cura_config_path(prof), overrides, cura_engine,
                                      log=log)
        else:
            log(f"[slice] using existing planar G-code {sliced_gcode}")
        retraction = 1.0
        if cura_info:
            retraction = float(cura_info["extruder"].get("retraction_amount",
                                                         cura_info["global"].get("retraction_amount")))
        else:  # sliced elsewhere: the mapper must know the file's retraction to recognise retract/unretract moves
            retraction = planar_retraction(planar_path, default=retraction)
            log(f"[slice] retraction length in that file: {retraction:g} mm")

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
            if mp.get("LIMIT_AXIS_SPEEDS"):
                from . import feed_limits
                with TIMER("axis speed limits"):
                    stats.update(feed_limits.apply(out_gcode, mp["MAX_SPEED_C"], mp["MAX_SPEED_B"], mp["MAX_SPEED_X"],
                                                   mp["MAX_SPEED_Z"], nozzle_offset=mp["NOZZLE_OFFSET"]))
        if support_check:
            with TIMER("4. support check"):
                from .params import expand_iterations
                support = sc.check(model_path, np.asarray(deformed.points), planar_path,
                                   part_offset=expand_iterations(p)[0], retraction_length=retraction,
                                   seg_size=mp["SEG_SIZE"])
            stats["floating_points"] = support["floating_points"]
            log(sc.format_report(support))
    poles = quality.vertical_extrusion(out_gcode, nozzle_offset=mp["NOZZLE_OFFSET"], retraction=retraction)
    stats["poles"] = len(poles[0])
    if support_check:
        stats["ungrounded_mm"] = round(support["ungrounded_mm"], 1)
        stats["ungrounded_real_mm"] = round(support["ungrounded_real_mm"], 1)
        stats["ungrounded_wall_mm"] = round(support["ungrounded_wall_mm"], 1)
        for k in ("island", "cantilever", "bridge"):  # floating runs by what their ends attach to
            stats[f"{k}_mm"] = round(support["runs"][k]["length_mm"], 1)
            stats[f"{k}_wall_mm"] = round(support["runs"][k]["wall_mm"], 1)  # of it walls/skin, not sparse infill
    # the time R-Theta Sim will show for this file (same planner; the settings header is comments only)
    from . import print_time
    with TIMER("print-time estimate"):
        try:
            stats["estimated_print_time_min"] = round(print_time.estimate(out_gcode, mp) / 60, 1)
        except (ValueError, TypeError) as e:  # a machine limit set to 0 or null in the profile
            log(f"[print-time] no estimate: {e}")

    # settings comments for the R-Theta simulator, at the top of the final G-code
    from . import sim_header
    from . import fast_deform as _fd
    info = cura_info
    if info is None:  # sliced elsewhere: rebuild the Cura settings (no slicing) so ;SETTING_3 can still be written
        try:
            from .cura import build_settings
            proj, g, e = build_settings(profiles.cura_config_path(prof), overrides, cura_engine)
            info = {"project": proj, "global": g, "extruder": e}
        except Exception:
            info = None
    shown = dict(prof); shown["map"] = dict(prof["map"], RETRACTION_LENGTH=retraction)  # the value actually used
    header = sim_header.build(name, shown, stats, planar_path, info, _fd.LAST_DEFORM_INFO,
                              sliced_by_cura_here=sliced_gcode is None)
    sim_header.prepend(out_gcode, header)
    log(quality.format_vertical(poles, line_offset=len(header)))  # line numbers in the final file
    total = time.perf_counter() - t0
    log(f"[map] {stats}")
    log(f"[done] {out_gcode}  ({total:.1f} s)")
    result = {"out": out_gcode, "stl": stl_path, "planar": planar_path, "stats": stats, "total_s": total,
              "timings": TIMER.as_dict()}
    with open(os.path.join(work_dir, f"timings_{impl}.json"), "w") as fh:
        json.dump(result, fh, indent=1, default=str)
    return result
