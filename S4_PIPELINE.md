# S4 headless pipeline

One command, STL in, 4-axis G-code out, for the Core R-Theta printer. The output uses C/X/Z/B/E moves, `G93`
inverse-time feed and `M83` relative E. It runs the same maths as `main.ipynb`, calls CuraEngine directly (no Cura
window), and is about 6× faster.

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
| `[deform]` | mesh size and where the deformed STL went |
| `[slice]` | the key Cura values actually used (layer height, line width, retraction, flavor, ...) |
| `[support]` | plastic that would be printed in mid-air, by kind (see below) |
| `[map]` | point counts from the mapping step |
| timing table | seconds per stage |

The **`[support]`** check looks for extruded points with nothing printed earlier (and not the bed) within 1 mm. It
sorts them into runs by what their ends connect to:

| kind | meaning | how bad |
|---|---|---|
| **island** | connected to nothing | truly unprintable, the one to watch |
| **cantilever** | anchored at one end | risky if long |
| **bridge** | anchored at both ends | normal FDM bridging |

Some bridging is normal. Sparse gyroid infill does it even in a flat print: a plain 30 mm cube shows 1.3%
"floating" infill with this check. Skip the check with `--no-support-check`.

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
| `ROTATION_ITERATIONS`, `DEFORMATION_ITERATIONS` | 100, 1000 | solver budgets. Changing them changes the result, not just the speed |
| `iterations` | (none) | multi-iteration schedule, see above |

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

**Compare runs using the same planar G-code** (`--sliced-gcode`). CuraEngine isn't deterministic: two slices of the
same STL differ by about 0.01 mm on a few hundred lines.

## 7. How it works

1. **Deform** (`s4/fast_deform.py`): tetgen, then the rotation-field and deformation least-squares solves.
2. **Slice** (`s4/cura.py`): headless CuraEngine.
3. **Map** (`s4/fast_map.py`): planar G-code mapped back to the 4-axis machine.
4. **Support check** (`s4/support_check.py`).

**Cura settings** are resolved from the 3mf the way the Cura GUI does it. Precedence is user > quality_changes >
quality > material > definition_changes > definition, and the extruder stack falls back to the global stack.
`=expressions`, `resolve` and `limit_to_extruder` are handled too. The result goes to CuraEngine as a generated
`.def.json`, because about 700 `-s` flags would exceed the Windows command-line limit.

**The fast path** gives exactly the reference's output. Its main changes, none of which alter the result:
- neighbour lists are rebuilt with numpy, but in VTK's exact order
- one Dijkstra instead of three identical ones
- one smoothing pass instead of 30 identical ones
- Jacobians built straight into CSR
- a row-parallel numba mat-vec inside scipy's TRF solver, with the same summation order
- batched LAPACK for the per-cell and per-point maths
- a regex G-code reader

## 8. Differences from the notebook

These are on by default. `--notebook-exact` turns all of them off and reproduces the notebook's output byte for byte.

| fix | notebook behaviour | now |
|---|---|---|
| **retraction "sticks"** (`SPLIT_RETRACTIONS`) | Each 1 mm unretract was a 1 mm plunge while extruding, and each retract a lift while retracting, at F20000. This left filament sticks at almost every travel (about 940 on the pi). | Retract in place, lift, travel, lower, unretract in place, at Cura's retraction speed (F3600). |
| **blotchy extrusion** (`SMOOTH_EXTRUSION_MULTIPLIER`, `EXTRUSION_MULTIPLIER_RANGE`) | The volume-compensation factor was constant per tetrahedron, so flow jumped at every tet boundary. On the pi there were 3,366 jumps of more than 25%. | Per-vertex factor blended smoothly along the path and clamped to 0.5×–2×: 84 jumps over 25%, total plastic within 0.2%. |
| **start-code blob** (`strip_start_prime`) | Cura's prime (`G1 F200 E3` at z = 20) was mapped into the part as a floating 3 mm blob. | Start code set to `G28 ; home`; the S4 header already primes at home. |
| **`NOZZLE_OFFSET`** | hard-coded 42 in the writer | taken from the profile |

## 9. Known limitations

- **Some plastic still prints in mid-air over flat bridges.** The deformation can't make a flat underside spanning a
  gap fully printable. The pi's crossbar is at z = 10 mm over an 11 mm gap, and its underside gets printed before the
  legs reach it. S4 is still far better than flat printing:

  | | pi | benchy |
  |---|---|---|
  | floating overall (flat print → S4) | 2.2% → 0.42% | 0.63% → 0.51% |
  | islands | 3 runs, 4 mm | 25 runs, 21 mm |

  Tuning the settings moved the pi between 0.18% and 0.77%, never to zero, and the best settings push B to its −130°
  limit. Fixing it properly means changing the rotation-field maths, which hasn't been done.
- **The deformation solve is quartic and unconverged.** Its residual `||N V − R N V0||²` gets squared again by least
  squares, and the solve stops at 1000 evaluations. That makes it chaotic: nudging the start point by 1e-9 mm moves
  vertices by up to 0.6 mm. So results can only be reproduced bit-for-bit, and they can differ slightly across
  machines or BLAS thread counts. An exact linear solve would take 0.2 s, but it would move benchy vertices by up to
  13 mm; not adopted.
- **The rotation smoothing is also quartic** (`W·Δ²` residuals, a `W²Δ⁴` penalty). Left unchanged.
- **`INITIAL_ROTATION_FIELD_SMOOTHING`** does one pass for any non-zero value. The notebook's loop recomputes from the
  same field every pass; it was probably meant to be iterative.
- **Old notebook G-code** (in `input_gcode/` / `output_gcode/` on older checkouts) was sliced in the Cura GUI with a different profile (0.1 mm layers,
  6.5 mm retraction). The old benchy used absolute extrusion, which the mapper can't handle.
