# S4 headless pipeline

One command, STL in, 4-axis G-code out, for the Core R-Theta printer. The output uses C/X/Z/B/E moves, `G93`
inverse-time feed and `M83` relative E. It's built on `main.ipynb`, calls CuraEngine directly (no Cura window), and
by default uses a deformation that avoids floating islands. What changed and why, with measurements, is in
[CHANGELOG.md](CHANGELOG.md).

**This page is the user guide.** The rest lives in separate documents:

| I want to... | Read |
|---|---|
| look up what a setting does | [docs/SETTINGS.md](docs/SETTINGS.md) |
| change the code, or check a change | [docs/DEVELOPING.md](docs/DEVELOPING.md) |
| view a print and its estimated time | [sim/README.md](sim/README.md) (R-Theta Sim) |
| see what changed, with measurements | [CHANGELOG.md](CHANGELOG.md) |

**Contents**
1. [Setup](#1-setup-once-per-computer) (once per computer)
2. [Slicing a model](#2-slicing-a-model)
3. [Tutorial: changing settings](#3-tutorial-changing-settings)
4. [Reading a run's output](#4-reading-a-runs-output)
5. [Known limitations](#5-known-limitations)

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

The pipeline writes only `output_gcode/` (results) and `build/` (work files). Both are git-ignored, so slicing never
shows up as a change to commit. Share G-code by sending the file itself, plus its `params_used.json` (section 3,
step 5).

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

What each setting does is in [docs/SETTINGS.md](docs/SETTINGS.md).

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
- **Compare settings with the `[support]` lines** (section 4). Less ungrounded plastic is better.
- **`--notebook-exact`** turns off every fix this fork adds and reproduces the notebook's output
  ([which fix is which setting](docs/SETTINGS.md#differences-from-the-notebook)).

## 4. Reading a run's output

A run prints these lines, in order:

| line | tells you |
|---|---|
| `[params]` | which profile file and overrides were used |
| `[deform]` | the deformation's own check: tilt achieved vs aimed for (`tilt_deg` / `target_tilt_deg`), folded tets and their volume, island seeds left. It prints a WARNING if the tilt falls well short |
| `[slice]` | the key Cura values actually used (layer height, line width, retraction, flavor, ...) |
| `[support]` | plastic that would be printed in mid-air, by kind (below) |
| `[quality]` | **poles**: extrusion dragged more than 2 mm straight down from a travel (should always say "none"); and extrusion running along the nozzle's own axis, which pushes into or pulls out of the bead. Sideways printing is normal S4 and isn't flagged |
| `[map]` | point counts from the mapping step |
| timing table | seconds per stage |

### The `[support]` lines
The check maps the toolpath into real space and walks it in print order. The first line is the one to watch:
**ungrounded** extrusion, meaning plastic with no chain of support down to the bed. Zero is the goal; a few tens of
mm of top skin or infill is normal (see [known limitations](#5-known-limitations)). The largest ungrounded regions
are listed with their layers and position, and the line also splits the length into walls/skin and sparse infill.
Judge by the walls/skin number, because sparse gyroid infill floats a little even in a flat print.

The lines after that are a finer, local view, sorted by what each floating run connects to:

| kind | meaning | how bad |
|---|---|---|
| **island** | connected to nothing | truly unprintable, the one to watch |
| **cantilever** | anchored at one end | risky if long |
| **bridge** | anchored at both ends | normal FDM bridging |

The exact rules are in [docs/DEVELOPING.md](docs/DEVELOPING.md#the-support-check). Skip the check with
`--no-support-check`.

### The settings block at the top of the G-code
Every output file starts with a block of comments: Cura's header lines, `;SETTING_3` lines (the Cura settings changed
from the profile defaults), and `; s4: KEY = value` lines (the model name and pipeline version, every deform and map
setting, the deformation's diagnostics, and the run's numbers such as `ungrounded_mm` and `poles`). It has no
timestamps, so two runs with the same settings give identical files. [R-Theta Sim](sim/README.md) reads it to show
and compare runs.

### Files in `build/<model>/`

| file | contents |
|---|---|
| `<model>_deformed_tet.stl` | the deformed mesh sent to Cura |
| `<model>_deformed_tet.gcode` | Cura's planar G-code |
| `cura_engine.log` | the exact CuraEngine command and its log |
| `params_used.json` | the complete settings for this run |
| `timings_fast.json` | stage timings and stats |
| `deformed_points.npy` | deformed vertex positions |

## 5. Known limitations

- **A little ungrounded plastic remains** (benchy: about 40 mm, mostly top skin and infill at the top of the hull;
  that's skin laid over sparse infill, which is normal bridging).
- **A part centred on the turntable** (all the test models are) needs fast C rotation where its paths pass near the
  axis: about 57°/r per mm of path. That is the machine's geometry, not a defect. The axis speed limits slow those
  moves down.
- **Squirtle has a 1.3 mm skin dot** at about (−10, −1, 31), in the shape the fit produces. Raising
  `ISLAND_LIFT_SLOPE` is stricter but distorts the part more.
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
