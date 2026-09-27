# Changelog

Changes in this fork relative to [jyjblrd/S4_Slicer](https://github.com/jyjblrd/S4_Slicer). Newest first. How to use the
pipeline is in [S4_PIPELINE.md](S4_PIPELINE.md). Measurements are on the `benchy upsidedown tilted` (33,156 tets) and
`pi 3mm` models with default settings, unless noted.

---

## 2026-09-27: Settings block for R-Theta Sim

### Added
- **Settings comments at the top of every output G-code** (`s4/sim_header.py`):
  - Cura's header lines
  - Cura-format `;SETTING_3` lines: every Cura setting changed from the profile defaults. They're generated for
    headless slices, which CuraEngine alone doesn't provide.
  - `; s4: KEY = value` lines: the settings used, the deformation's diagnostics, and the run's quality numbers

  The block is deterministic, so fast and reference outputs stay identical. `--notebook-exact` output matches the
  notebook byte for byte after the block.

### Notes
- `;SETTING_3` is written exactly as Cura writes it: the JSON is escaped (`\` → `\`, newline → `
`) before
  it's split into 80-character comment lines. A reader has to undo that escaping before parsing the JSON, as Cura's
  own G-code reader does.
- CuraEngine's placeholder `;MINX…;MAXZ` values (`2.14748e+06`) and `;Filament used: 0m` aren't copied, since they
  aren't real measurements.

---

## 2026-09-27: No more "poles" at the base

### Fixed (on by default; `--notebook-exact` restores the notebook's behaviour)
- **Vertical filament "poles"** (`SAFE_TRAVEL_TRANSITIONS`). When a travel leaves the part's mesh, the notebook lifts
  it to the highest point printed so far, and keeps the first point back on the part at that height too. If the next
  move prints, it extruded from the lifted point straight down to the layer. On the upside-down benchy's first
  layers, that drew 8 sticks of 4–6 mm standing on the bed. Now the nozzle travels down to the re-entry point's true
  height first. If the re-entry move is itself a print move, it's travelled over instead of extruded, because its
  start was off the part.
- **Extruding `G00` lines.** When a move needs a big B rotation, it's split into small steps. The notebook labelled
  those steps with the *previous* move's command, so the first print move after a travel came out as `G00` lines
  carrying extrusion. Each step now carries its own move's command.

### Added
- **`[quality]` line after every run.** It rebuilds the nozzle-tip path from the final 4-axis G-code and reports
  **poles**: extruding moves that go more than 2 mm straight down right after a travel. It also reports steep
  extruding segments for information.

### Results
- **Benchy (same planar slice):** poles 8 → 0; steep extruding segments 7 → 1 (the one left is a real steep stretch
  of layer).
- **Every model in `input_models/`** (benchy, pi, Squirtle, dino, z mount, B-axis mount) reports no poles.

---

## 2026-09-27: Island-free deformation

### Added
- **`DEFORMATION_METHOD: "island_free"`, the new default** (`s4/island_free.py`). It replaces the step that turns the
  rotation (tilt) field into a deformed shape; the rotation field itself is still the notebook's.
  - **Fold-free fit:** each tet is fitted to its target rotation, plus a barrier term that becomes infinite before
    any tet can turn inside out. Every step is capped before the first tet would flip, and the tilt is applied in
    10 stages, because large rotations can't be reached in one step. The solver is L-BFGS, preconditioned with the
    mesh stiffness matrix, with numba kernels.
  - **Island lifting:** a priority-flood from the bed finds every vertex that can't be reached from the bed by a path
    rising at least `ISLAND_LIFT_SLOPE` per mm. Those vertices get soft height targets, and the fit is solved again.
  - New settings: `ISLAND_LIFT_SLOPE`, `ISLAND_LIFT_ROUNDS`, `LIFT_WEIGHT`, `BARRIER_WEIGHT`, `FLIP_FREE_STAGES`,
    `FLIP_FREE_STAGE_ITERATIONS`. `--notebook-exact` selects the old method.
- **Grounded support report.** The `[support]` report now leads with *ungrounded* extrusion: plastic with no chain of
  support down to the bed. The old check only counted an island's first layer, because each later layer looked
  "supported" by the floating one below it. `support_check.check_planar()` runs the same check on Cura's planar
  toolpath.

### Results

| | notebook method | `island_free` |
|---|---|---|
| benchy: ungrounded extrusion in Cura's toolpath | 1,258 mm (two islands of 624 and 239 mm) | 29 mm (largest: 17 mm of top skin) |
| pi: ungrounded extrusion in Cura's toolpath | 18 mm | 0 mm |
| benchy: inverted tets | 202 | 0 |
| benchy: extrusion points needing more than 2× plastic | 3.8% | 0.4% |
| benchy: deformation time | ~60 s | ~40 s |
| benchy, 5-iteration recipe: ungrounded extrusion | – | ~51 mm |

### Investigation
- **Where the islands come from.** Cura slices the deformed shape flat, so every local low point of that shape (lower
  than its surroundings, not on the bed) starts printing in mid-air. Such low points, measured by how long they float
  before joining the bed-connected part, match the toolpath islands exactly. The benchy's two islands are its two
  largest low points.
- **Why the notebook leaves them.** Its deformation solve stops at 1000 evaluations, unconverged. It misses its own
  target rotations by 9° median (26° at the 95th percentile), and folds 202 tets. Near the islands the target tilt is
  also too small: flat downward faces need about 60°, and the smoothing settles at about 43°.
- **Joshua's original has them too.** His committed pi toolpath has 17 mm of ungrounded extrusion, about the same as
  the notebook method here. His vertex order couldn't be matched to ours (a different tetgen version), so this was
  measured on his planar toolpath directly.
- **Tried and rejected:**
  - parameter tuning: the pi moved between 0.18% and 0.77% floating, never to zero, and B hit its −130° limit
  - moving the part off-centre: no help
  - an exact linear deformation solve: better fit, but more floating extrusion (2,000 mm)
  - quadratic rotation smoothing: worse
  - per-cell overhang feedback: worse
  - lifting heights directly: inverted thousands of tets
  - pinning the bottom face flat: jammed the solver

---

## 2026-09-26: Output fixes, build profiles, docs

### Added
- **Build profiles** (`s4/profile.py`): every setting for a build in one JSON file, with `deform`, `map` and `cura`
  sections. `params/<model name>.json` is picked up automatically.
  - Command-line options: `--set KEY=VALUE`, `--show-params`, `--init-params`.
  - Misspelled keys are rejected with a suggestion; angle keys can be given in degrees (`_DEG`).
  - Every run saves `params_used.json`, which reproduces the run byte for byte.
  - `params/pi 3mm.json` is the default profile for the pi.
- **Support check** (`s4/support_check.py`): mid-air extrusion reported after every run, sorted into island,
  cantilever and bridge runs.
- **Settings tutorial** in `S4_PIPELINE.md`, plus `requirements.txt` with the verified package versions.

### Fixed (on by default; `--notebook-exact` restores the notebook's behaviour)
- **Filament "sticks" at every retraction** (`SPLIT_RETRACTIONS`). The notebook treats points between a 1 mm retract
  and the matching unretract as travel, and lifts travel 1 mm along the tool axis. So each unretract was written as a
  1 mm plunge while extruding 1 mm of filament, about 100× a normal segment, at F20000. With the 3mf's 1 mm
  retraction this happened about 940 times on the pi. Old GUI slices used 6.5 mm, which never matched the notebook's
  1.0 constant. The pipeline now retracts in place, lifts, travels, lowers, and unretracts in place, at Cura's
  retraction speed (F3600).
- **Blotchy extrusion** (`SMOOTH_EXTRUSION_MULTIPLIER`, `EXTRUSION_MULTIPLIER_RANGE`). The volume-compensation factor
  was constant per tet, so flow jumped at tet boundaries: 3,366 jumps of more than 25% on the pi. It's now a
  per-vertex factor, blended along the path and clamped to 0.5×–2×. That leaves 84 such jumps, and total plastic
  changes by 0.2%.
- **Floating blob from Cura's start code** (`strip_start_prime`). The prime (`G1 F200 E3` at z = 20) was mapped into
  the part as a 3 mm blob in mid-air. The start code is now `G28 ; home`.
- **`NOZZLE_OFFSET`** was hard-coded to 42 in the G-code writer; it now comes from the profile.
- **Misspelled Cura setting names** used to be silently ignored by CuraEngine; they're now rejected before the run.

### Changed
- **G-code flavor: RepRap** (was Marlin).
- **`main.ipynb`:** rotation limits of ±360° (the committed notebook still had ±3600°), GIFs off, a diagnostic cell
  for cells printing in air, and a commented-out tetgen `maxvolume` experiment.
- **Generated outputs are no longer tracked:** `gifs/`, `input_gcode/`, `output_gcode/`, `output_models/` and
  `pickle_files/` keep only a `.gitkeep`.
- **Docs** reorganized into `S4_PIPELINE.md`, with a pointer from the README.

### Removed
- Dead code: `CSRPattern`, the old `deg:` parameter hack, `tools/test_fast_map.py`.

### Investigation
- **The floating plastic over the pi's bridge** comes from the slicer, not the simulator. Cura's own planar toolpath
  already floats in the deformed shape.
- **Sparse gyroid infill "floats" by nature:** a plain 30 mm cube shows 1.3% under a 1 mm-radius check. Bridges
  anchored at both ends are normal.
- **CuraEngine isn't deterministic.** Two slices of the same STL differ by about 0.01 mm on a few hundred lines, even
  single-threaded. Compare runs on the same planar G-code.

---

## 2026-09-25: Headless pipeline

### Added
- **`s4_slice.py`:** STL in, 4-axis G-code out, with no Cura window.
- **Headless CuraEngine** (`s4/cura.py`). It resolves the `cura_config.3mf` setting stack the way the Cura GUI does:
  container precedence, `=expressions`, `resolve` and `limit_to_extruder`. The result goes to CuraEngine as a
  generated `.def.json`, because about 700 `-s` flags would exceed the Windows command-line limit. The enforced
  settings are origin at plate centre, no Z-hop, relative extrusion and a centred model.
- **Reference implementation** (`s4/reference.py`): a line-for-line port of `main.ipynb`. It's verified
  **byte-identical** to the notebook's own benchy G-code, and its deformed mesh is **bit-identical** to the notebook's
  pickle.
- **Fast implementation** (`s4/fast_deform.py`, `s4/fast_map.py`, `s4/fast_lsq.py`), byte-identical to the
  reference:
  - neighbour lists rebuilt with numpy, but in VTK's exact order
  - one Dijkstra instead of three identical ones
  - one smoothing pass instead of 30 identical ones
  - Jacobians built straight into CSR
  - a row-parallel numba mat-vec inside scipy's TRF solver, with the same summation order
  - batched LAPACK for the per-cell and per-point maths
  - a regex G-code reader instead of pygcode
  - skipping the unused `shared_vertices` loop and the recomputation of deformed attributes
- **Equivalence gate** (`tools/check_equivalence.py`) and 4-axis G-code comparison (`tools/compare_gcode.py`).
- **Multi-iteration schedules** replaying the notebook's cell 4 → 7 → 9 loop; `params/benchy_upsidedown_tilted_recipe.json`
  holds the recipe from notebook cell 5.

### Performance (benchy, laptop on battery)

| stage | notebook | fast |
|---|---|---|
| mesh setup | 57 s | 5 s |
| rotation field | 81 s | 11 s |
| deformation solve | 165 s | ~60 s |
| mapping + writing | 131 s | 6 s |
| end to end, same run | 858 s | 131 s |

### Investigation
- **The notebook's deformation solve is quartic and unconverged,** and chaotic: nudging the start point by 1e-9 mm moves
  vertices by up to 0.6 mm. So it can only be reproduced bit-for-bit, and a GPU port couldn't match it within a small
  tolerance. It was sped up on the CPU instead.
- **The rotation smoothing is also quartic** (`W·Δ²` residuals). Left unchanged, and still unchanged today.
- **`INITIAL_ROTATION_FIELD_SMOOTHING = 30` does one pass:** each of the 30 recomputes from the same unsmoothed field.
- **The old GUI-sliced benchy used absolute extrusion (M82).** The mapper treats E as relative, so that output had
  per-segment values like `E1794`.
