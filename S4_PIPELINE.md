# S4 headless pipeline

One command, STL in, 4-axis G-code out, for the Core R-Theta printer. The output uses C/X/Z/B/E moves, `G93`
inverse-time feed and `M83` relative E. It's built on `main.ipynb`, calls CuraEngine directly (no Cura window), and
by default uses a deformation that avoids floating islands. What changed and why, with measurements, is in
[CHANGELOG.md](CHANGELOG.md).

**Contents**
1. [Setup](#1-setup-once-per-computer) (once per computer)
2. [Slicing a model](#2-slicing-a-model)
3. [Tutorial: changing settings](#3-tutorial-changing-settings)
4. [Reading a run's output](#4-reading-a-runs-output)
5. [Settings reference](#5-settings-reference)
6. [Checking code changes](#6-checking-code-changes)
7. [How it works](#7-how-it-works)
8. [Differences from the notebook](#8-differences-from-the-notebook)
9. [Known limitations](#9-known-limitations)

All commands run from the `S4_Slicer` folder, in PowerShell or cmd.

---

## 1. Setup (once per computer)

1. Install **Python 3.14** and **UltiMaker Cura 5.13**. The pipeline finds `CuraEngine.exe` under
   `C:\Program Files\UltiMaker Cura *`; otherwise set the `CURA_ENGINE` environment variable or pass `--cura-engine`.
2. Create the virtual environment and install the packages:
   ```
   py -3.14 -m venv venv
   venv\Scripts\python -m pip install -r requirements.txt
   ```
3. Check that it works. This prints the settings without slicing anything:
   ```
   venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --show-params
   ```

## 2. Slicing a model

```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl"
```
This deforms the model, slices it with CuraEngine, maps it back to 4 axes, and writes
`output_gcode/<model>.gcode`. Like the notebook, it **overwrites** an existing file with that name; use
`-o other.gcode` to choose another. A pi takes about 15 s and a benchy about 2 min.

Generated files stay out of git: the contents of `output_gcode/`, `input_gcode/`, `output_models/`, `gifs/`,
`pickle_files/`, `build/` and `artifacts/` are ignored, so slicing never shows up as a change to commit.
Share G-code by sending the file itself, plus its `params_used.json` (section 3, step 5).

```
venv\Scripts\python s4_slice.py "input_models/benchy upsidedown tilted.stl" -o my_benchy.gcode
venv\Scripts\python s4_slice.py "input_models/benchy upsidedown tilted.stl" --params params/benchy_upsidedown_tilted_recipe.json
venv\Scripts\python s4_slice.py --help
```
The second line runs the benchy's multi-iteration recipe from notebook cell 5.

## 3. Tutorial: changing settings

There are three ways, from quickest to most permanent:

| I want to... | Do this |
|---|---|
| try a value once | add `--set KEY=VALUE` (deform/map settings) or `--cura-set KEY=VALUE` (Cura settings) |
| keep settings for a model | put them in `params/<model name>.json`; it's used automatically every time |
| keep several variants of one model | put each in its own file and pick one with `--params FILE` |

### Step 1: see the current settings
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --show-params
```
The first line says where the settings came from: `built-in defaults`, or the profile file that was picked up.
Nothing is sliced.

### Step 2: try a change once
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --set MAX_OVERHANG=10 --cura-set layer_height=0.1
```
Repeat `--set` / `--cura-set` as often as you like. Values are read as JSON, so lists work:
`--set "PART_OFFSET=[0,10,0]"`. Nothing is saved; the next run is back to normal.

### Step 3: save the settings for this model
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --init-params
```
This writes `params/pi 3mm.json` with every setting at its current value (plus any `--set` you added). It refuses to
overwrite an existing file unless you add `--force`. Open the file, change what you want, and **delete the lines you
didn't change**, so it shows only what's special about this build:
```json
{
  "description": "pi: gentler tilt, 0.1 mm layers",
  "deform": {"MAX_OVERHANG": 20, "NEIGHBOUR_LOSS_WEIGHT": 50},
  "cura":   {"set": {"layer_height": 0.1, "infill_sparse_density": 15}}
}
```
From now on the plain command uses it, and the first output line confirms that:
```
[params] ...\params\pi 3mm.json (matched model name)
```
Commit the file so the rest of the team slices the model the same way.

### Step 4: keep variants side by side
Copy the profile to a new name and choose it explicitly. Use `-o` so the variants don't overwrite each other:
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --params "params/pi 3mm - fine.json" -o "output_gcode/pi fine.gcode"
```

### Step 5: repeat or share an exact run
Every run saves the complete settings it used to `build/<model>/params_used.json`. Pass that file back to
reproduce the run exactly, or send it along with a G-code file:
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --params "build/pi 3mm/params_used.json"
```

### Good to know
- **Precedence**, later wins: defaults < profile file < `--set` < `--cura-set` / `--cura-config` / `--notebook-exact`.
- **Typos are caught before anything runs**, with a suggestion:
  `error: unknown deform setting 'MAX_OVERHNG' (did you mean MAX_OVERHANG?)`. This works for Cura names too.
- **`PART_OFFSET` is subtracted**, as in the notebook: `[0, 10, 0]` moves the part 10 mm toward −Y.
- **Compare settings with the `[support]` lines** (see section 4). Fewer islands is better.

## 4. Reading a run's output

A run prints these lines, in order:

| line | tells you |
|---|---|
| `[params]` | which profile file and overrides were used |
| `[deform]` | the deformation's own check: tilt achieved vs aimed for (`tilt_deg` / `target_tilt_deg`), folded tets and their volume, island seeds left. It prints a WARNING if the tilt falls well short |
| `[slice]` | the key Cura values actually used (layer height, line width, retraction, flavor, ...) |
| `[support]` | plastic that would be printed in mid-air, by kind (see below) |
| `[quality]` | poles: extrusion dragged more than 2 mm straight down from a travel (should always say "none"), and extrusion running along the nozzle's own axis (pushing into or pulling out of the bead). Vertical moves with the nozzle tilted sideways are normal S4 printing and aren't flagged |
| `[map]` | point counts from the mapping step |
| timing table | seconds per stage |

The **`[support]`** check maps the toolpath into real space and walks it in print order. Its first line is the one
to watch: **ungrounded** extrusion, meaning plastic with no chain of support down to the bed. A point counts as
grounded if the bed or grounded plastic from an earlier layer is within 1 mm, or if it's within 5 mm along the same
extrusion line of such a point (normal bridging and overhang reach). Everything else is ungrounded: a floating
island, plus every layer stacked on top of it. The largest ungrounded regions are listed with their layers and
position.

The lines after that are a finer, local view. They count points with nothing printed earlier within 1 mm, sorted by
what the ends of each run connect to:

| kind | meaning | how bad |
|---|---|---|
| **island** | connected to nothing | truly unprintable, the one to watch |
| **cantilever** | anchored at one end | risky if long |
| **bridge** | anchored at both ends | normal FDM bridging |

Some bridging is normal: sparse gyroid infill does it even in a flat print. Skip the check with `--no-support-check`.

**The top of every output G-code** is a settings block, which R-Theta Sim reads to show and compare runs:
- **Cura's header lines:** `;FLAVOR:`, `;Layer height:`, `;Generated with …`.
- **`;SETTING_3` lines** in Cura's own format: every Cura setting changed from the profile defaults. They're copied
  from the planar G-code when it came from the Cura GUI, and generated for headless slices, because CuraEngine alone
  doesn't write them.
- **`; s4: KEY = value` lines:**
  - the model name and pipeline version
  - every deform and map setting (per iteration, too), and the Cura overrides
  - the deformation's diagnostics
  - the run's numbers: `failed_points`, `floating_points`, `ungrounded_mm`, `poles`, ...

The block is deterministic (no timestamps), so two runs with the same settings give identical files.

Files in `build/<model>/`:

| file | contents |
|---|---|
| `<model>_deformed_tet.stl` | the deformed mesh sent to Cura |
| `<model>_deformed_tet.gcode` | Cura's planar G-code |
| `cura_engine.log` | the exact CuraEngine command and its log |
| `params_used.json` | the complete settings for this run |
| `timings_fast.json` | stage timings and stats |
| `deformed_points.npy` | deformed vertex positions |

## 5. Settings reference

### Profile file format
Every section and key is optional; anything left out keeps its default.
```json
{
  "description": "free text",
  "deform": {"MAX_OVERHANG": 10, "PART_OFFSET": [0, 0, 0],
             "iterations": [{"NEIGHBOUR_LOSS_WEIGHT": 100}, {"NEIGHBOUR_LOSS_WEIGHT": 50}]},
  "map":    {"NOZZLE_OFFSET": 41.5, "MIN_ROTATION": -120},
  "cura":   {"config": "cura_config.3mf", "set": {"layer_height": 0.1}, "strip_start_prime": true}
}
```
- **Iterations:** `iterations` replays the notebook's "run cell 4 → 7 → 9 again" loop. Each entry lists only what
  changes, and inherits everything else from the previous iteration.
- **Degrees:** the rotation-limit keys the notebook kept in radians can be written in degrees with a `_DEG` suffix:
  `MAX_POS_ROTATION_DEG`, `MAX_NEG_ROTATION_DEG`, `ROTATION_MAX_DELTA_DEG`.
- **Old files:** a flat dict, or `{"iterations": [...]}`, still loads and is read as the `deform` section.

### `deform`: how the part is warped (notebook cells 2–7)
| setting | default | what it does |
|---|---|---|
| `MAX_OVERHANG` | 30 | overhang angle (°) the deformation aims for. Lower = warps harder, more tilt |
| `ROTATION_MULTIPLIER` | 2 | scales the target tilt. Higher = more tilt |
| `NEIGHBOUR_LOSS_WEIGHT` | 30 | smoothness of the tilt field. Higher = smoother, gentler changes |
| `SET_INITIAL_ROTATION_TO_ZERO` | false | pull areas without overhangs toward no tilt (less noisy) |
| `STEEP_OVERHANG_COMPENSATION` | true | extra rotation for cells that would print in air |
| `INITIAL_ROTATION_FIELD_SMOOTHING` | 30 | 0 = off; any other value = one smoothing pass (see section 9) |
| `MAX_POS_ROTATION`, `MAX_NEG_ROTATION` | ±360° | clamp on the target rotation |
| `PART_OFFSET` | [0, 0, 0] | shift the part on the plate, in mm (subtracted) |
| `ROTATION_ITERATIONS`, `DEFORMATION_ITERATIONS` | 100, 1000 | solver budgets (`DEFORMATION_ITERATIONS` only applies to `notebook`). Changing them changes the result, not just the speed |
| `iterations` | (none) | multi-iteration schedule, see above |
| `DEFORMATION_METHOD` | `island_free` | how the tilt field becomes a deformed shape: `island_free` (section 7) or `notebook` |
| `ISLAND_LIFT_SLOPE` | 0.5 | `island_free`: every point must be reachable from the bed rising at least this much per mm. Higher = stricter (1.0 ≈ 45° overhangs) but more distortion |
| `ISLAND_LIFT_ROUNDS` | 5 | `island_free`: rounds of lifting |
| `FIT_METHOD` | `penalty` | `island_free`: `penalty` (robust on any mesh) or `barrier` (strictly fold-free, but can stall on fine meshes) |
| `FOLD_PENALTY`, `FOLD_PENALTY_DET` | 100, 0.5 | `penalty`: how hard tets squashed below 0.5× volume (or folded) are pushed back. Lower lets thin features squash (the pi then extrudes along the nozzle axis) |
| `PENALTY_ITERATIONS`, `LIFT_ITERATIONS` | 300, 100 | `penalty`: solver iterations for the fit, and per lifting round |
| `BED_PIN_WEIGHT`, `BED_TOL` | 0, 0.3 | `island_free`: the bed face (vertices within `BED_TOL` mm of the lowest point) always counts as supported. `BED_PIN_WEIGHT` > 0 also holds it flat on the bed; off by default, since it costs tilt and support (benchy: 1475 mm unsupported at 500) |
| `LIFT_WEIGHT` | 50 | `island_free`: strength of the lift targets. Too weak and the fold penalty wins over the lift, leaving islands (Squirtle kept an unsupported tower at 5) |
| `FLIP_FREE_STAGES`, `FLIP_FREE_STAGE_ITERATIONS`, `BARRIER_WEIGHT` | 10, 150, 0.02 | `barrier` only: tilt ramp stages and barrier strength |
| `SLIVER_QUALITY`, `MICRO_TET_VOLUME` | 0.03, 1e-3 | `barrier` only: badly shaped or tiny tets get no barrier (they would stall it) |

### `map`: the machine and the final 4-axis moves (notebook cells 17–18)
| setting | default | what it does |
|---|---|---|
| `NOZZLE_OFFSET` | 42 | mm from the B pivot to the nozzle tip (the notebook notes 41.5 is the true value) |
| `MIN_ROTATION`, `MAX_ROTATION` | −130, 30 | B-axis limits in degrees |
| `ROTATION_AVERAGING_ALPHA` | 0.2 | smoothing of B along the path. Lower = smoother, slower to react |
| `ROTATION_MAX_DELTA` | 1° | B steps bigger than this are split into smaller moves |
| `SEG_SIZE` | 0.6 | mm; long moves are split into segments this long |
| `MAX_EXTRUSION_MULTIPLIER` | 10 | hard cap on the extrusion compensation |
| `EXTRUSION_MULTIPLIER_RANGE` | [0.5, 2.0] | clamp on the extrusion compensation (`null` = off) |
| `SMOOTH_EXTRUSION_MULTIPLIER` | true | smooth the extrusion compensation (section 8) |
| `SPLIT_RETRACTIONS` | true | retract/unretract in place (section 8) |
| `SAFE_TRAVEL_TRANSITIONS` | true | after a travel that left the part, lower before printing; rotation steps keep their own move's command (section 8) |
| `LIMIT_AXIS_SPEEDS` | true | slow down any move that would drive an axis past the limits below (section 8) |
| `MAX_SPEED_C`, `MAX_SPEED_B` | 360, 180 | deg/s, bed rotation and nozzle tilt. Placeholders: set your machine's real limits |
| `MAX_SPEED_X`, `MAX_SPEED_Z` | 150, 50 | mm/s. Placeholders, as above |
| `RETRACTION_LENGTH` | null | mm; `null` = use Cura's `retraction_amount` |

### `cura`: slicing
| setting | default | what it does |
|---|---|---|
| `config` | `cura_config.3mf` | Cura project with the base profile (save one from Cura to use your own) |
| `set` | {} | Cura setting overrides by **internal** name, e.g. `layer_height`, `infill_sparse_density`, `infill_pattern`, `wall_line_count`, `speed_print`, `material_print_temperature` |
| `strip_start_prime` | true | drop Cura's start-code prime (section 8) |

The full list of Cura names is in
`C:\Program Files\UltiMaker Cura 5.13.0\share\cura\resources\definitions\fdmprinter.def.json`.
Whatever the profile says, these are always enforced, because the mapper depends on them:

| setting | value | why |
|---|---|---|
| machine_center_is_zero | True | origin at the build-plate centre |
| retraction_hop_enabled | False | no Z-hop |
| machine_gcode_flavor | RepRap (RepRap) | the printer's firmware flavor |
| relative_extrusion | True | the mapper scales each move's E value |
| center_object | True | model centred on the plate ("Arrange All") |

## 6. Checking code changes

After any change to `s4/`, run the equivalence gate:
```
venv\Scripts\python tools\check_equivalence.py "input_models/pi 3mm.stl"
venv\Scripts\python tools\check_equivalence.py "input_models/pi 3mm.stl" --params params/test_two_iterations.json
```
It runs the reference implementation (a line-for-line port of the notebook, `--impl reference`) and the fast one on
the same model and the same planar G-code. It compares the deformed mesh and the final G-code, and exits non-zero on
failure. Today they are byte-identical. Add `--notebook-exact` to check the notebook's original behaviour.

Other tools:
```
venv\Scripts\python tools\compare_gcode.py a.gcode b.gcode
venv\Scripts\python s4_slice.py model.stl --notebook-exact
venv\Scripts\python tools\run_reference.py map --model M.stl --out DIR --sliced planar.gcode --deformed-pkl pickle_files/deformed_M.pkl
```
- `compare_gcode.py` compares any two 4-axis files, with tolerances.
- `--notebook-exact` reproduces old notebook results.
- `run_reference.py` maps a notebook pickle.

**Compare runs using the same planar G-code** (`--sliced-gcode`), because CuraEngine isn't deterministic. After a
change that affects output, add an entry to [CHANGELOG.md](CHANGELOG.md).

## 7. How it works

1. **Deform** (`s4/fast_deform.py`, `s4/island_free.py`): tetgen, the notebook's rotation (tilt) field, then turning
   that field into a deformed shape.
2. **Slice** (`s4/cura.py`): headless CuraEngine.
3. **Map** (`s4/fast_map.py`): planar G-code mapped back to the 4-axis machine.
4. **Support check** (`s4/support_check.py`).

**The deformation (`island_free`, default).** Cura slices the deformed shape flat. Any local low point of that shape,
a spot lower than everything around it that isn't on the bed, starts printing in mid-air. `island_free` keeps the
notebook's tilt field but builds the shape in two steps:
- **fit:** first the exact least-squares fit of every tet to its target rotation (a single sparse solve), then a soft
  penalty pushes open any tet that was squashed or folded (`FIT_METHOD: penalty`). `FIT_METHOD: barrier` forbids
  folds outright instead, but it can stall on fine meshes.
- **lifting:** vertices that can't be reached from the bed by a path rising at least `ISLAND_LIFT_SLOPE` per mm get
  height targets, and the fit is solved again. Those regions then print later, growing out from where they're
  attached.

`DEFORMATION_METHOD: notebook` uses the notebook's own least-squares solve instead. Before/after numbers are in
the [changelog](CHANGELOG.md#2026-09-27-island-free-deformation).

**Cura settings** are resolved from the 3mf the way the Cura GUI does it. Precedence is user > quality_changes >
quality > material > definition_changes > definition, and the extruder stack falls back to the global stack.
`=expressions`, `resolve` and `limit_to_extruder` are handled too. The result goes to CuraEngine as a generated
`.def.json`, because about 700 `-s` flags would exceed the Windows command-line limit.

**Two implementations.** `--impl fast` (default) and `--impl reference`, a line-for-line port of the notebook, give
byte-identical output (section 6). The speed-ups in the fast one are listed in the
[changelog](CHANGELOG.md#2026-09-25-headless-pipeline).

## 8. Differences from the notebook

These are on by default. `--notebook-exact` turns all of them off and reproduces the notebook's output byte for byte
(everything after the settings block described in section 4).
The rotation (tilt) field is the notebook's in every case.

| what | setting | difference |
|---|---|---|
| deformation | `DEFORMATION_METHOD` | fold-free, island-free shape instead of the notebook's unconverged solve |
| retractions | `SPLIT_RETRACTIONS` | retract/unretract in place, at Cura's retraction speed (the notebook extruded during a 1 mm plunge, leaving "sticks") |
| extrusion compensation | `SMOOTH_EXTRUSION_MULTIPLIER`, `EXTRUSION_MULTIPLIER_RANGE` | blended smoothly along the path and clamped to 0.5×–2× (was constant per tet, so flow jumped) |
| start code | `strip_start_prime` | Cura's prime is dropped (it became a floating blob inside the part) |
| travel re-entry | `SAFE_TRAVEL_TRANSITIONS` | after a travel that left the part, the nozzle lowers before printing (the notebook printed downward from the lifted point, drawing "poles"); rotation-split steps keep their own move's command; points mapped just below the bed are clamped to it instead of dropped |
| nozzle offset | `NOZZLE_OFFSET` | taken from the profile (was hard-coded to 42) |
| axis speeds | `LIMIT_AXIS_SPEEDS`, `MAX_SPEED_*` | each move's G93 time was its planar move's time, which ignores what the axes do: a 1° tilt swings X/Z ~0.7 mm, moves near the centre need big C rotations, travels lifted over the part are much longer than planned. Up to 900,000°/s on C. Now each move takes at least as long as every axis needs at its limit; the G94 `F20000` moves become G93 too |

Why each of these was needed, with measurements, is in the [changelog](CHANGELOG.md).

## 9. Known limitations

- **A little ungrounded plastic remains** (benchy: about 30–50 mm, mostly top skin and infill at the top of the hull;
  that's skin laid over sparse infill, which is normal bridging). Raising `ISLAND_LIFT_SLOPE` is stricter but distorts
  the part more.
- **`island_free` changes the shape Cura slices**, so its output differs from the notebook's (by design). Vertex
  positions are deterministic run to run.
- **The notebook's deformation solve** (`DEFORMATION_METHOD: notebook`) is unconverged and chaotic. Its results can
  only be reproduced bit-for-bit, and can differ slightly across machines or BLAS thread counts.
- **The rotation smoothing is quartic** (`W·Δ²` residuals, a `W²Δ⁴` penalty). Left unchanged.
- **`INITIAL_ROTATION_FIELD_SMOOTHING`** does one pass for any non-zero value. The notebook's loop recomputes from the
  same field every pass; it was probably meant to be iterative.
- **Cura GUI slices need relative extrusion** (M83). The mapper treats every E value as relative, so a slice made with
  absolute extrusion (M82) comes out badly over-extruded. Headless slicing always uses relative.
- **CuraEngine isn't deterministic:** compare runs on the same planar G-code (`--sliced-gcode`).
