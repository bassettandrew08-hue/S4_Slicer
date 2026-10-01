"""
Settings comments at the top of the final 4-axis G-code, for the R-Theta simulator (and for people):

1. Cura's header lines (;FLAVOR, ;Layer height, ;Generated with, ...), copied from the planar G-code.
2. ;SETTING_3 lines in Cura's own format (a JSON object, chunked into comment lines): every Cura setting changed
   from the profile defaults. Copied from the planar G-code when Cura wrote them (GUI slices); generated from the
   3mf's user/quality_changes containers plus the pipeline's overrides for headless slices, because CuraEngine on
   its own doesn't write them.
3. "; s4: KEY = value" lines: the S4 settings used and the run's quality numbers.

Everything here is deterministic (no timestamps, no implementation name), so the fast and reference pipelines
still produce identical files.
"""
import json
import math
import os
import re
import subprocess

from .profile import DEG_KEYS  # angle settings, shown as deg2rad(<degrees>)

HEAD = re.compile(r"^;\s*(FLAVOR|TIME|Filament used|Layer height|MINX|MINY|MINZ|MAXX|MAXY|MAXZ|Generated with|"
                  r"TARGET_MACHINE\.NAME|PRINT\.TIME|NOZZLE_DIAMETER|MATERIAL)(\s*[:=]|\s)", re.I)
SENTINEL = "2.14748e+06"   # CuraEngine's placeholder for MINX..MAXZ when not run from the Cura GUI
PREFIX = ";SETTING_3 "


def cura_header_lines(planar_path):
    """Cura's header comments from the top of the planar G-code (before the first command)."""
    out = []
    with open(planar_path, errors="replace") as fh:
        for line in fh:
            t = line.strip()
            if not t:
                continue
            if not t.startswith(";"):
                break
            if HEAD.match(t) and SENTINEL not in t and t != ";Filament used: 0m":
                out.append(t)
    return out


def setting3_from_file(planar_path):
    with open(planar_path, errors="replace") as fh:
        return [l.rstrip("\n") for l in fh if l.startswith(PREFIX)]


def _inst_cfg(name, definition, values, extra_meta=None):
    lines = ["[general]", "version = 4", f"name = {name}", f"definition = {definition}", "",
             "[metadata]", "type = quality_changes", "quality_type = fast", "setting_version = 23"]
    for k, v in (extra_meta or {}).items():
        lines.append(f"{k} = {v}")
    lines += ["", "[values]"]
    for k in sorted(values):
        v = values[k]
        if isinstance(v, bool):
            v = "True" if v else "False"
        v = str(v).replace("\n", "\n\t")  # multi-line values (start g-code) as in Cura's .cfg files
        lines.append(f"{k} = {v}")
    return "\n".join(lines) + "\n\n"


def setting3_generated(cura_info):
    """Cura-format ;SETTING_3 lines for a headless slice (the settings changed from the definition defaults)."""
    proj = cura_info["project"]
    g_vals, e_vals = {}, {}
    # container order is highest priority first: apply lowest first so higher ones win
    for cont, target in ((proj.gstack.containers, g_vals), (proj.estack.containers, e_vals)):
        user, quality_changes = (cont + [{}, {}])[:2]
        target.update(quality_changes)
        target.update(user)
    for k, v in proj.overrides_global.items():
        if not proj.settable_per_extruder(k):
            g_vals[k] = v
    for k, v in proj.overrides_extruder.items():
        if proj.settable_per_extruder(k):
            e_vals[k] = v
    data = {"global_quality": _inst_cfg("s4 pipeline", "custom", g_vals),
            "extruder_quality": [_inst_cfg("s4 pipeline", "custom", e_vals, {"position": 0})]}
    s = json.dumps(data)
    # Cura escapes these before splitting into comment lines (GCodeWriter._serialiseSettings)
    s = s.replace("\\", "\\\\").replace("\n", "\\n").replace("\r", "\\r")
    width = 80 - len(PREFIX)
    return [PREFIX + s[i:i + width] for i in range(0, len(s), width)]


def _fmt(v):
    if isinstance(v, str):
        return v.replace("\r", "\\r").replace("\n", "\\n")  # a raw newline would end the comment line
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return f"{v:g}"
    if isinstance(v, (list, tuple)):
        return ", ".join(_fmt(x) for x in v)
    return str(v)


def _fmt_param(k, v):
    if k in DEG_KEYS and isinstance(v, (int, float)):
        return f"deg2rad({math.degrees(v):g})"
    return _fmt(v)


def git_version():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        r = subprocess.run(["git", "describe", "--always", "--dirty"], cwd=here, capture_output=True, text=True,
                           timeout=5)
        return r.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def s4_lines(model_name, profile, stats, deform_info=None, planar_source=None):
    """The '; s4: KEY = value' lines: model, pipeline version, every deform setting (and per-iteration overrides),
    every map setting, the Cura config and overrides, the deformation diagnostics (without timings, so the file
    stays deterministic) and the run's stats. Angles are shown as deg2rad(<degrees>).
    """
    from .params import expand_iterations
    out = [f"; s4: model = {model_name}", f"; s4: pipeline_version = {git_version()}"]
    if profile.get("description"):
        out.append(f"; s4: description = {_fmt(profile['description'])}")
    if planar_source:
        out.append(f"; s4: planar_gcode = {planar_source}")
    deform = profile["deform"]
    _, iters = expand_iterations(deform)
    out.append(f"; s4: iterations = {len(iters)}")
    for k, v in deform.items():
        if k == "iterations":
            continue
        out.append(f"; s4: {k} = {_fmt_param(k, v)}")
    if len(iters) > 1:
        for i, it in enumerate(deform["iterations"], 1):
            for k, v in it.items():
                out.append(f"; s4: iteration_{i}.{k} = {_fmt_param(k, v)}")
    out.append("; s4: tetrahedralize = default (no maxvolume)")
    for k, v in profile["map"].items():
        out.append(f"; s4: {k} = {_fmt_param(k, v)}")
    out.append(f"; s4: cura_config = {profile['cura']['config']}")
    for k, v in sorted(profile["cura"]["set"].items()):
        out.append(f"; s4: cura_set.{k} = {_fmt(v)}")
    if deform_info:
        for k, v in deform_info.items():
            if k in ("seconds", "fit_seconds"):
                continue  # keep the header deterministic
            out.append(f"; s4: deform.{'solver_iterations' if k == 'iterations' else k} = {_fmt(v)}")
    stat_names = {"lost_vertices": "failed_points", "gcode_points": "planar_points",
                  "bad_barycentric_sum": "bad_barycentric_points"}
    for k, v in stats.items():
        out.append(f"; s4: {stat_names.get(k, k)} = {_fmt(v)}")
    return out


def build(model_name, profile, stats, planar_path, cura_info=None, deform_info=None, sliced_by_cura_here=True):
    lines = ["; ---- slicer settings (read by R-Theta Sim) ----"]
    lines += cura_header_lines(planar_path)
    s3 = setting3_from_file(planar_path)
    if not s3 and cura_info is not None:
        s3 = setting3_generated(cura_info)
    lines += s3
    lines += s4_lines(model_name, profile, stats, deform_info,
                      planar_source=None if sliced_by_cura_here else os.path.basename(planar_path))
    lines.append("; ---- end of slicer settings ----")
    return lines


def prepend(path, lines):
    with open(path, "r") as fh:
        body = fh.read()
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n" + body)
