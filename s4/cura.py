"""
Headless CuraEngine slicing driven by a Cura project file (cura_config.3mf).

Cura's GUI resolves every setting (definition defaults, '=expression' values, and
the user/quality/material containers) before handing them to CuraEngine. CuraEngine
itself only knows `default_value`s, so we do that resolution here:

  global stack   : user > quality_changes > intent > quality > material > variant
                   > definition_changes > machine definition (custom -> fdmprinter)
  extruder stack : same order for the extruder's own containers, then falls back to
                   the global stack. Settings with settable_per_extruder=False are
                   always read from the global stack. '=expressions' are evaluated in
                   the context of the stack that asked (as Cura's SettingFunction does).

The resolved global values are written into a generated .def.json (passed with -j,
which sidesteps the 32k-character Windows command-line limit); extruder values that
differ from the global ones are passed with -e0 -s key=value.

Reference for the approach: github.com/austinsus/universal-s4-slicer (GPL-3.0).
"""
import ast
import configparser
import glob
import json
import math
import os
import re
import subprocess
import tempfile
import xml.etree.ElementTree as ET
import zipfile

CURA_SEARCH = [
    r"C:\Program Files\UltiMaker Cura *\CuraEngine.exe",
    r"C:\Program Files\Ultimaker Cura *\CuraEngine.exe",
]

# Settings the S4 inverse mapper relies on; these are enforced after resolution.
S4_REQUIRED = {
    "machine_center_is_zero": True,      # origin at build-plate centre
    "retraction_hop_enabled": False,     # no Z-hop
    "machine_gcode_flavor": "RepRap (RepRap)",  # RepRap flavour (what the printer firmware uses)
    "relative_extrusion": True,          # mapper scales per-move E values
    "center_object": True,               # model centred on the plate ("Arrange All")
}

# XML material setting key -> Cura setting key (subset of Cura's XmlMaterialProfile map)
MATERIAL_XML_MAP = {
    "print temperature": "default_material_print_temperature",
    "heated bed temperature": "default_material_bed_temperature",
    "standby temperature": "material_standby_temperature",
    "processing temperature graph": "material_flow_temp_graph",
    "print cooling": "cool_fan_speed",
    "retraction amount": "retraction_amount",
    "retraction speed": "retraction_speed",
    "adhesion tendency": "material_adhesion_tendency",
    "surface energy": "material_surface_energy",
    "build volume temperature": "build_volume_temperature",
    "anti ooze retract position": "material_anti_ooze_retracted_position",
    "anti ooze retract speed": "material_anti_ooze_retraction_speed",
    "break preparation position": "material_break_preparation_retracted_position",
    "break preparation speed": "material_break_preparation_speed",
    "break preparation temperature": "material_break_preparation_temperature",
    "break position": "material_break_retracted_position",
    "flush purge speed": "material_flush_purge_speed",
    "flush purge length": "material_flush_purge_length",
    "end of filament purge speed": "material_end_of_filament_purge_speed",
    "end of filament purge length": "material_end_of_filament_purge_length",
    "maximum park duration": "material_maximum_park_duration",
    "no load move factor": "material_no_load_move_factor",
    "break speed": "material_break_speed",
    "break temperature": "material_break_temperature",
    "shrinkage percentage": "material_shrinkage_percentage",
}

_MISSING = object()


def find_cura_engine():
    env = os.environ.get("CURA_ENGINE")
    if env and os.path.exists(env):
        return env
    hits = []
    for pattern in CURA_SEARCH:
        hits += glob.glob(pattern)
    if not hits:
        raise FileNotFoundError("CuraEngine.exe not found; pass --cura-engine or set CURA_ENGINE")
    return sorted(hits)[-1]


class Definition:
    """Flattened setting definitions: name -> property dict."""

    def __init__(self):
        self.props = {}

    def _walk(self, node):
        for name, spec in node.items():
            if not isinstance(spec, dict):
                continue
            if spec.get("type") != "category":
                self.props.setdefault(name, {}).update({k: v for k, v in spec.items() if k != "children"})
            if "children" in spec:
                self._walk(spec["children"])

    def load_chain(self, path, search_dirs):
        d = json.load(open(path, encoding="utf-8"))
        if "inherits" in d:
            parent = _find_def(d["inherits"], search_dirs)
            self.load_chain(parent, search_dirs)
        self._walk(d.get("settings", {}))
        for name, spec in d.get("overrides", {}).items():
            self.props.setdefault(name, {}).update(spec)
        return self

    def get(self, name, prop, default=None):
        return self.props.get(name, {}).get(prop, default)


def _find_def(name, dirs):
    for d in dirs:
        p = os.path.join(d, name + ".def.json")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"definition {name}.def.json not found in {dirs}")


def _parse_typed(text, typ):
    """Convert an .inst.cfg literal string to the setting's Python type."""
    if not isinstance(text, str):
        return text
    s = text.strip()
    if typ == "bool":
        return s.lower() in ("true", "1", "yes", "on")
    if typ == "int" or typ == "extruder" or typ == "optional_extruder":
        try:
            return int(float(s))
        except ValueError:
            return s
    if typ == "float":
        try:
            return float(s.replace(",", "."))
        except ValueError:
            return s
    if typ in ("[int]", "polygon", "polygons"):
        try:
            return ast.literal_eval(s)
        except (ValueError, SyntaxError):
            return s
    return text  # str / enum: keep verbatim (start g-code keeps its newlines)


class Stack:
    def __init__(self, name, containers, definition, resolver, next_stack=None):
        self.name = name
        self.containers = containers  # list of dicts, highest priority first
        self.definition = definition  # Definition for this stack's bottom container
        self.resolver = resolver
        self.next_stack = next_stack
        self.cache = {}


class CuraProject:
    CONTAINER_ORDER = ["user", "quality_changes", "intent", "quality", "material", "variant", "definition_changes"]

    def __init__(self, threemf_path, cura_engine=None):
        self.engine = cura_engine or find_cura_engine()
        res = os.path.join(os.path.dirname(self.engine), "share", "cura", "resources")
        self.def_dirs = [os.path.join(res, "definitions"), os.path.join(res, "extruders")]
        self.overrides_global = {}
        self.overrides_extruder = {}
        self._resolving = set()
        self._resolving_resolve = set()
        self.failed = set()
        self._load_3mf(threemf_path)

    # ---------------------------------------------------------------- loading
    def _load_3mf(self, path):
        with zipfile.ZipFile(path) as zf:
            names = zf.namelist()
            cfgs = {}
            for n in names:
                if n.startswith("Cura/") and (n.endswith(".cfg")):
                    cp = configparser.ConfigParser(interpolation=None)
                    cp.optionxform = str
                    cp.read_string(zf.read(n).decode("utf-8"))
                    cfgs[n] = cp
            materials = {n: zf.read(n).decode("utf-8") for n in names if n.endswith(".xml.fdm_material")}

        # instance containers by id (file name without extension; also their [general] name)
        instances = {}
        global_cfg = extruder_cfg = None
        for n, cp in cfgs.items():
            base = os.path.basename(n)
            if n.endswith(".global.cfg"):
                global_cfg = cp
            elif n.endswith(".extruder.cfg"):
                extruder_cfg = cp
            elif n.endswith(".inst.cfg"):
                instances[base[:-len(".inst.cfg")]] = cp
        if global_cfg is None or extruder_cfg is None:
            raise ValueError(f"{path}: no global/extruder stack found (is this a Cura project file?)")

        material_values = {}
        for n, xml in materials.items():
            material_values[os.path.basename(n)[:-len(".xml.fdm_material")]] = _parse_material_xml(xml)

        def build(cp):
            ids = [cp["containers"][k] for k in sorted(cp["containers"], key=int)]
            definition_id = ids[-1]
            containers = []
            for cid in ids[:-1]:
                if cid in instances and instances[cid].has_section("values"):
                    containers.append(dict(instances[cid]["values"]))
                elif cid in material_values:
                    containers.append(material_values[cid])
                else:
                    containers.append({})  # empty_* containers
            return definition_id, containers

        gdef_id, gcont = build(global_cfg)
        edef_id, econt = build(extruder_cfg)
        self.global_def = Definition().load_chain(_find_def(gdef_id, self.def_dirs), self.def_dirs)
        self.extruder_def = Definition().load_chain(_find_def(edef_id, self.def_dirs), self.def_dirs)
        self.gstack = Stack("global", gcont, self.global_def, self)
        self.estack = Stack("extruder", econt, self.extruder_def, self, next_stack=self.gstack)
        self.machine_name = global_cfg["general"].get("name", "?")

    # ------------------------------------------------------------- resolution
    def set(self, key, value, extruder=True, global_=True):
        """Force a value (like a user-container change) on the global and/or extruder stack."""
        if global_:
            self.overrides_global[key] = value
        if extruder:
            self.overrides_extruder[key] = value
        self.gstack.cache.clear()
        self.estack.cache.clear()

    def prop(self, key, name, default=None):
        v = self.extruder_def.get(key, name, _MISSING)
        if v is _MISSING:
            v = self.global_def.get(key, name, default)
        return v

    def settable_per_extruder(self, key):
        v = self.prop(key, "settable_per_extruder", True)
        return v if isinstance(v, bool) else str(v).lower() == "true"

    def _type(self, key):
        return self.prop(key, "type", "str")

    def _raw(self, stack, key):
        """First container (or definition) value in stack order; returns (raw, is_expression)."""
        overrides = self.overrides_extruder if stack is self.estack else self.overrides_global
        if key in overrides:
            return overrides[key], False
        for c in stack.containers:
            if key in c:
                raw = c[key]
                if isinstance(raw, str) and raw.startswith("="):
                    return raw[1:], True
                return _parse_typed(raw, self._type(key)), False
        if key in stack.definition.props:
            spec = stack.definition.props[key]
            if "value" in spec:
                v = spec["value"]
                if isinstance(v, str):
                    return v, True
                return v, False
            if "default_value" in spec:
                return spec["default_value"], False
        if stack.next_stack is not None:
            return self._raw(stack.next_stack, key)
        return None, False

    def value(self, stack, key):
        if key in stack.cache:
            return stack.cache[key]
        if stack is self.estack and not self.settable_per_extruder(key):
            v = self.value(self.gstack, key)
            stack.cache[key] = v
            return v
        guard = (stack.name, key)
        if guard in self._resolving:
            return self.prop(key, "default_value")
        self._resolving.add(guard)
        try:
            if stack is self.gstack:
                # GlobalStack: 'resolve' property wins unless the user container sets the key.
                # While a resolve is being evaluated, nested lookups of the same key skip it
                # and fall through to the normal container lookup (Cura's _resolving_settings).
                resolve = self.global_def.get(key, "resolve")
                user_has = key in self.overrides_global or key in (stack.containers[0] if stack.containers else {})
                if resolve and not user_has and key not in self._resolving_resolve:
                    self._resolving_resolve.add(key)
                    self._resolving.discard(guard)
                    try:
                        v = self._eval(stack, key, resolve)
                    finally:
                        self._resolving_resolve.discard(key)
                    stack.cache[key] = v
                    return v
                # GlobalStack: settings with a limit_to_extruder (even -1 -> active extruder)
                # are read from that extruder; single-extruder machine -> extruder 0
                if self.global_def.get(key, "limit_to_extruder") is not None and self.settable_per_extruder(key):
                    v = self.value(self.estack, key)
                    stack.cache[key] = v
                    return v
            raw, is_expr = self._raw(stack, key)
            v = self._eval(stack, key, raw) if is_expr else raw
        finally:
            self._resolving.discard(guard)
        stack.cache[key] = v
        return v

    def _eval(self, stack, key, expr):
        proj = self

        def extruder_value(_pos, k):
            return proj.value(proj.estack, k)

        env = {
            "math": math,
            "re": re,
            "extruderValue": extruder_value,
            "extruderValues": lambda k: [proj.value(proj.estack, k)],
            "anyExtruderWithMaterial": lambda k: 0,
            "anyExtruderNrWithOrDefault": lambda k: 0,
            "resolveOrValue": lambda k: proj.value(proj.gstack, k),
            "defaultExtruderPosition": lambda: "0",
            "valueFromContainer": lambda k, i: proj.value(stack, k),
            "valueFromExtruderContainer": lambda k, i: proj.value(proj.estack, k),
            "extruderPosition": lambda *a: 0,
        }
        try:
            tree = ast.parse(expr.strip(), mode="eval")
            for node in ast.walk(tree):
                if (isinstance(node, ast.Name) and node.id not in env
                        and (node.id in self.global_def.props or node.id in self.extruder_def.props)):
                    env[node.id] = self.value(stack, node.id)
            return eval(compile(tree, key, "eval"), env)  # noqa: S307 (trusted Cura profile)
        except Exception:
            self.failed.add(key)
            return self.prop(key, "default_value")

    def resolve_all(self):
        keys = sorted(set(self.global_def.props) | set(self.extruder_def.props))
        g = {}
        e = {}
        for k in keys:
            gv = self.value(self.gstack, k) if k in self.global_def.props else None
            if gv is not None:
                g[k] = gv
            if self.settable_per_extruder(k):
                ev = self.value(self.estack, k)
                if ev is not None:
                    e[k] = ev
        # values StartSliceJob computes on the fly
        start = re.sub(r";.+?(\n|$)", "\n", str(g.get("machine_start_gcode", "")))
        bed = ["material_bed_temperature", "material_bed_temperature_layer_0"]
        g["material_bed_temp_prepend"] = re.search(r"\{(%s)(,\s?\w+)?\}" % "|".join(bed), start) is None
        prt = ["material_print_temperature", "material_print_temperature_layer_0", "default_material_print_temperature",
               "material_initial_print_temperature", "material_final_print_temperature", "material_standby_temperature",
               "print_temperature"]
        g["material_print_temp_prepend"] = re.search(r"\{(%s)(,\s?\w+)?\}" % "|".join(prt), start) is None
        for code_key in ("machine_start_gcode", "machine_end_gcode"):
            if code_key in g:
                g[code_key] = _expand_tokens(str(g[code_key]), {**g, **e})
        return g, e


def _expand_tokens(text, values):
    def rep(m):
        k = m.group(1)
        return str(values[k]) if k in values else m.group(0)
    return re.sub(r"\{(\w+)\}", rep, text)


def _parse_material_xml(xml_text):
    """Top-level (machine-independent) settings of an .xml.fdm_material file."""
    root = ET.fromstring(xml_text)
    ns = {"m": root.tag.split("}")[0].strip("{")} if root.tag.startswith("{") else {}
    q = (lambda t: f"m:{t}") if ns else (lambda t: t)
    out = {}
    settings = root.find(q("settings"), ns)
    if settings is not None:
        for s in settings.findall(q("setting"), ns):  # direct children only (not per-machine)
            key = MATERIAL_XML_MAP.get(s.get("key"))
            if key:
                out[key] = s.text.strip() if s.text else ""
    props = root.find(q("properties"), ns)
    if props is not None:
        dia = props.find(q("diameter"), ns)
        if dia is not None and dia.text:
            out["material_diameter"] = dia.text.strip()
    return out


def _fmt(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, dict, tuple)):
        return json.dumps(v)
    return str(v)


def check_setting_names(proj_or_3mf, names, cura_engine=None):
    """Raise ValueError for Cura setting names that don't exist (CuraEngine would silently ignore them)."""
    import difflib
    proj = proj_or_3mf if isinstance(proj_or_3mf, CuraProject) else CuraProject(proj_or_3mf, cura_engine)
    known = set(proj.global_def.props) | set(proj.extruder_def.props)
    for k in (names or {}):
        if k not in known:
            m = difflib.get_close_matches(k, known, n=1)
            raise ValueError(f"unknown Cura setting {k!r}" + (f" (did you mean {m[0]}?)" if m else "")
                             + ". Use Cura's internal setting name, e.g. layer_height, infill_sparse_density")


def build_settings(threemf, overrides=None, cura_engine=None):
    """Resolve the project and apply S4-required + user overrides. Returns (project, global, extruder)."""
    proj = CuraProject(threemf, cura_engine)
    check_setting_names(proj, overrides)
    for k, v in {**S4_REQUIRED, **(overrides or {})}.items():
        if isinstance(v, str) and k in proj.global_def.props:
            v = _parse_typed(v, proj._type(k)) if proj._type(k) not in ("str", "enum") else v
        proj.set(k, v)
    g, e = proj.resolve_all()
    return proj, g, e


def slice_stl(stl_path, out_gcode, threemf, overrides=None, cura_engine=None, log=print, threads=None):
    proj, g, e = build_settings(threemf, overrides, cura_engine)
    for k, want in S4_REQUIRED.items():
        got = e.get(k, g.get(k))
        if got != want:
            raise RuntimeError(f"setting {k} resolved to {got!r}, S4 needs {want!r}")

    tmp = tempfile.mkdtemp(prefix="s4_cura_")
    def_path = os.path.join(tmp, "s4_resolved.def.json")
    children = {}
    for k, v in g.items():
        children[k] = {"label": k, "type": "str", "default_value": _fmt(v)}
    json.dump({"version": 2, "name": "s4_resolved", "metadata": {},
               "settings": {"s4": {"label": "s4", "type": "category", "children": children}}},
              open(def_path, "w", encoding="utf-8"))

    args = [proj.engine, "slice", "-v", "-j", def_path]
    if threads:
        args.append(f"-m{threads}")
    args.append("-e0")
    ext_diff = {k: v for k, v in e.items() if _fmt(g.get(k, _MISSING)) != _fmt(v)}
    for k, v in ext_diff.items():
        args += ["-s", f"{k}={_fmt(v)}"]
    args += ["-l", os.path.abspath(stl_path), "-o", os.path.abspath(out_gcode)]
    cmdlen = sum(len(a) + 3 for a in args)
    if cmdlen > 30000:
        raise RuntimeError(f"CuraEngine command line too long ({cmdlen} chars)")

    key_info = {k: e.get(k, g.get(k)) for k in (
        "layer_height", "layer_height_0", "line_width", "wall_line_count", "infill_sparse_density",
        "infill_pattern", "retraction_amount", "retraction_hop_enabled", "relative_extrusion",
        "machine_gcode_flavor", "machine_center_is_zero", "center_object", "adhesion_type",
        "support_enable", "material_print_temperature", "speed_print", "speed_wall_0")}
    log(f"[slice] {proj.machine_name}: {len(g)} global + {len(ext_diff)} extruder-specific settings; "
        + ", ".join(f"{k}={v}" for k, v in key_info.items()))
    if proj.failed:
        log(f"[slice] note: {len(proj.failed)} expressions fell back to defaults: {', '.join(sorted(proj.failed)[:10])}")

    proc = subprocess.run(args, capture_output=True, text=True, errors="replace")
    errors = [ln for ln in proc.stderr.splitlines() if "[error]" in ln.lower()]
    if proc.returncode != 0 or not os.path.exists(out_gcode) or os.path.getsize(out_gcode) == 0:
        raise RuntimeError(f"CuraEngine failed (exit {proc.returncode}):\n" + proc.stderr[-4000:])
    if errors:
        log(f"[slice] CuraEngine reported {len(errors)} error lines, first: {errors[0]}")
    with open(os.path.join(os.path.dirname(os.path.abspath(out_gcode)), "cura_engine.log"), "w") as fh:
        fh.write(" ".join(f'"{a}"' if " " in a else a for a in args) + "\n\n" + proc.stderr)
    return {"global": g, "extruder": e, "args": args}
