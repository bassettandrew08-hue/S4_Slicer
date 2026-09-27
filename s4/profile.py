"""
Build profiles: every setting for one build in one JSON file.

    {
      "description": "pi at 0.1 mm layers",
      "deform": {"MAX_OVERHANG": 30, "PART_OFFSET": [0, 0, 0],
                 "iterations": [{"NEIGHBOUR_LOSS_WEIGHT": 100}, {"NEIGHBOUR_LOSS_WEIGHT": 50}]},
      "map":    {"NOZZLE_OFFSET": 42, "MIN_ROTATION": -130},
      "cura":   {"config": "cura_config.3mf", "set": {"layer_height": 0.1}}
    }

Every section and key is optional; anything left out keeps its default. Resolution order, later wins:
    built-in defaults  <  profile file  <  --set KEY=VALUE  <  --cura-set / --cura-config / --notebook-exact

The profile file is --params PATH, or params/<model name>.json if that file exists. Angle keys may also be
given in degrees by adding _DEG (MAX_POS_ROTATION_DEG, ROTATION_MAX_DELTA_DEG, ...). Older params files
(a flat dict of deformation keys, or {"iterations": [...]}) still work and are read as the "deform" section.
"""
import copy
import difflib
import json
import math
import os

from .params import DEFAULT_PARAMS

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# notebook cell 17/18 constants + the output fixes (see S4_PIPELINE.md)
MAP_DEFAULTS = dict(
    SEG_SIZE=0.6,                   # mm, planar moves are split into segments of at most this length
    MAX_ROTATION=30,                # deg, B-axis limit (positive tilt)
    MIN_ROTATION=-130,              # deg, B-axis limit (negative tilt)
    NOZZLE_OFFSET=42,               # mm, B pivot to nozzle tip
    ROTATION_AVERAGING_ALPHA=0.2,   # exponential smoothing of B along the path
    ROTATION_MAX_DELTA=math.radians(1),  # rad; B steps bigger than this get interpolated
    MAX_EXTRUSION_MULTIPLIER=10,
    RETRACTION_LENGTH=None,         # mm; None = take Cura's resolved retraction_amount
    SPLIT_RETRACTIONS=True,         # retract/unretract in place instead of during the travel lift/plunge
    SMOOTH_EXTRUSION_MULTIPLIER=True,
    EXTRUSION_MULTIPLIER_RANGE=[0.5, 2.0],  # or null for no clamp
)
CURA_DEFAULTS = dict(
    config="cura_config.3mf",       # Cura project with the base profile (relative to the repo folder)
    set={},                         # Cura setting overrides, e.g. {"layer_height": 0.1}
    strip_start_prime=True,         # drop Cura's start-code prime (it would land inside the part)
)
NOTEBOOK_EXACT = {"deform": {"DEFORMATION_METHOD": "notebook"},
                  "map": {"SPLIT_RETRACTIONS": False, "SMOOTH_EXTRUSION_MULTIPLIER": False,
                          "EXTRUSION_MULTIPLIER_RANGE": None},
                  "cura": {"strip_start_prime": False}}
DEG_KEYS = {"MAX_POS_ROTATION", "MAX_NEG_ROTATION", "ROTATION_MAX_DELTA"}
SECTIONS = ("deform", "map", "cura")


class ProfileError(ValueError):
    pass


def defaults():
    return {"description": "", "deform": dict(DEFAULT_PARAMS), "map": dict(MAP_DEFAULTS),
            "cura": copy.deepcopy(CURA_DEFAULTS)}


def _suggest(key, known):
    m = difflib.get_close_matches(key, known, n=1)
    return f" (did you mean {m[0]}?)" if m else ""


def _check_keys(section, d, known, where):
    for k in d:
        if k not in known:
            raise ProfileError(f"{where}: unknown {section} setting {k!r}{_suggest(k, known)}")


def _normalise_deg(d):
    """Convert KEY_DEG entries to radians under KEY."""
    out = {}
    for k, v in d.items():
        if k.endswith("_DEG") and k[:-4] in DEG_KEYS:
            out[k[:-4]] = math.radians(v)
        else:
            out[k] = v
    return out


def merge(profile, overlay, where="profile"):
    """Overlay a (possibly partial) profile dict onto a full one, validating keys."""
    if not isinstance(overlay, dict):
        raise ProfileError(f"{where}: must be a JSON object")
    overlay = {k: v for k, v in overlay.items() if not k.startswith("_")}  # "_comment" keys are ignored
    if not any(k in overlay for k in SECTIONS + ("description",)):
        overlay = {"deform": overlay}  # legacy flat params file
    for k in overlay:
        if k not in SECTIONS + ("description",):
            raise ProfileError(f"{where}: unknown section {k!r}{_suggest(k, SECTIONS)}")
    p = copy.deepcopy(profile)
    if "description" in overlay:
        p["description"] = overlay["description"]
    if "deform" in overlay:
        d = _normalise_deg({k: v for k, v in overlay["deform"].items() if not k.startswith("_")})
        its = d.pop("iterations", None)
        _check_keys("deform", d, list(DEFAULT_PARAMS), where)
        p["deform"].update(d)
        if its is not None:
            if not isinstance(its, list) or not its:
                raise ProfileError(f"{where}: deform.iterations must be a non-empty list of objects")
            clean = []
            for i, it in enumerate(its):
                it = _normalise_deg({k: v for k, v in it.items() if not k.startswith("_")})
                _check_keys("deform", it, [k for k in DEFAULT_PARAMS if k != "PART_OFFSET"], f"{where} iteration {i + 1}")
                clean.append(it)
            p["deform"]["iterations"] = clean
    if "map" in overlay:
        m = _normalise_deg({k: v for k, v in overlay["map"].items() if not k.startswith("_")})
        _check_keys("map", m, list(MAP_DEFAULTS), where)
        p["map"].update(m)
    if "cura" in overlay:
        c = {k: v for k, v in overlay["cura"].items() if not k.startswith("_")}
        _check_keys("cura", c, list(CURA_DEFAULTS), where)
        if "set" in c:
            p["cura"]["set"] = {**p["cura"]["set"], **c.pop("set")}
        p["cura"].update(c)
    return p


def _parse_value(text):
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        low = text.strip().lower()
        if low in ("true", "false"):
            return low == "true"
        return text


def apply_set(profile, assignment):
    """--set KEY=VALUE. KEY is a deform or map setting (optionally prefixed deform./map./cura.), VALUE is JSON."""
    key, sep, val = assignment.partition("=")
    if not sep:
        raise ProfileError(f"--set expects KEY=VALUE, got {assignment!r}")
    key, value = key.strip(), _parse_value(val.strip())
    section, _, bare = key.rpartition(".")
    bare_base = bare[:-4] if bare.endswith("_DEG") else bare
    if not section:
        if bare_base in DEFAULT_PARAMS:
            section = "deform"
        elif bare_base in MAP_DEFAULTS:
            section = "map"
        elif bare in CURA_DEFAULTS:
            section = "cura"
        else:
            raise ProfileError(f"--set: unknown setting {bare!r}{_suggest(bare, list(DEFAULT_PARAMS) + list(MAP_DEFAULTS))}"
                               " (Cura settings go in --cura-set)")
    return merge(profile, {section: {bare: value}}, where="--set")


def default_profile_path(model_path):
    name = os.path.splitext(os.path.basename(model_path))[0]
    return os.path.join(HERE, "params", f"{name}.json")


def resolve(model_path, params_path=None, sets=(), cura_sets=None, cura_config=None, notebook_exact=False):
    """Build the effective profile. Returns (profile, source_description)."""
    p = defaults()
    source = "built-in defaults"
    path = params_path or default_profile_path(model_path)
    if params_path or os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            try:
                data = json.load(fh)
            except json.JSONDecodeError as e:
                raise ProfileError(f"{path}: invalid JSON ({e})") from None
        p = merge(p, data, where=os.path.basename(path))
        source = path if params_path else f"{path} (matched model name)"
    for s in sets:
        p = apply_set(p, s)
    if cura_sets:
        p["cura"]["set"] = {**p["cura"]["set"], **cura_sets}
    if cura_config:
        p["cura"]["config"] = cura_config
    if notebook_exact:
        p = merge(p, NOTEBOOK_EXACT, where="--notebook-exact")
    return p, source


def deform_params(p):
    """The dict s4.params.expand_iterations / both deform implementations take."""
    return copy.deepcopy(p["deform"])


def cura_config_path(p):
    c = p["cura"]["config"]
    return c if os.path.isabs(c) else os.path.join(HERE, c)


def to_json(p):
    """Human-editable JSON; rotation limits also shown in degrees for readability."""
    out = copy.deepcopy(p)
    for sec in ("deform", "map"):
        for k in list(out[sec]):
            if k in DEG_KEYS:
                out[sec][k + "_DEG"] = round(math.degrees(out[sec].pop(k)), 6)
        for it in out[sec].get("iterations", []) if sec == "deform" else []:
            for k in list(it):
                if k in DEG_KEYS:
                    it[k + "_DEG"] = round(math.degrees(it.pop(k)), 6)
    text = json.dumps(out, indent=2)
    # keep short lists of numbers on one line ("PART_OFFSET": [0, 10, 0])
    import re
    return re.sub(r"\[\s*([-\d.eE+,\s]+?)\s*\]", lambda m: "[" + ", ".join(x.strip() for x in m.group(1).split(",")) + "]", text)
