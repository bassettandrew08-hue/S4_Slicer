# Developing the pipeline

For people changing the code. Using the pipeline is in [S4_PIPELINE.md](../S4_PIPELINE.md); every setting is in
[SETTINGS.md](SETTINGS.md). All commands run from the `S4_Slicer` folder.

**Contents**
- [The rules](#the-rules)
- [Checking code changes](#checking-code-changes)
- [How it works](#how-it-works)
- [Where the time goes](#where-the-time-goes)
- [Modules](#modules)

## The rules

1. **Fast and reference must be byte-identical.** `--impl fast` (default) and `--impl reference` (a line-for-line
   port of the notebook) must produce the same deformed mesh and the same G-code, byte for byte. Run the equivalence
   gate after any change to `s4/`.
2. **`--notebook-exact` must keep reproducing the notebook.** It switches off every fix this fork adds
   ([the list](SETTINGS.md#differences-from-the-notebook)); the output after the settings block must match the
   notebook's byte for byte. A new fix needs a setting that `--notebook-exact` turns off.
3. **Byrd's math stays untouched.** The rotation (tilt) field and its smoothing are Joshua Bird's, and are copied
   as they are, including the quartic smoothing. Changes replace the step *after* the field (the deformed shape) or
   the steps around the mapping, never the field itself. The fast code may only be faster, never different.
   One exception, opt-in and off by default: `PIT_FREE_ROUNDS` (below) changes the field where the fit leaves pits.
4. **Run the six-model suite before changing a default.** A fix for one model must not hurt the others.
5. **Compare runs on the same planar G-code** (`--sliced-gcode`), because CuraEngine isn't deterministic.
6. **Log what you changed.** After a change that affects output, add an entry to [CHANGELOG.md](../CHANGELOG.md),
   with the measurements that justify it.

## Checking code changes

### The equivalence gate
```
venv\Scripts\python tools\check_equivalence.py "input_models/pi 3mm.stl"
venv\Scripts\python tools\check_equivalence.py "input_models/pi 3mm.stl" --params params/test_two_iterations.json
```
It runs the reference implementation and the fast one on the same model and the same planar G-code, compares the
deformed mesh and the final G-code, and exits non-zero on failure. The two must be byte-identical (apart from the
`planar_gcode` line of the settings block). Add `--notebook-exact` to check the notebook's original behaviour, and
`--set KEY=VALUE` (repeatable, passed to both runs) to check an opt-in setting, e.g.
`--set MAP_SMOOTH=0.5 --set MAX_C_STEP=3`. `--extra ...` (last) goes to the candidate run only.

### The six-model suite
`tools/model_suite.py` slices all six test models (in parallel) and reports structure (ungrounded, island,
cantilever: walls/skin and total), poles, extrusion along the nozzle axis, tilt, roughness (the sim's path-roughness
count, `s4/quality.py`), the shape metrics of `tools/shape_metrics.py` (`outside_mm`, `surface_rough_pct`,
`asym_mm`), `waviness` (joints turning > 10 deg where Cura's planar path runs straight: what the mapping adds,
`tools/roughness_breakdown.py`) and print time. With `--baseline DIR` (an earlier `--out`) it flags regressions.
`--sliced-from DIR` reuses an earlier run's planar G-code, for comparing map-only settings without CuraEngine's
noise. `tools/roughness_breakdown.py RUN_DIR MODEL` sorts a run's rough segments by cause.
```
venv\Scripts\python tools\model_suite.py --out build/suite/base
venv\Scripts\python tools\model_suite.py --out build/suite/try --set KEY=VALUE --baseline build/suite/base
```
The roughness count depends on point spacing, so compare it only between runs with the same `SEG_SIZE`.

### Other tools
| tool | what it does |
|---|---|
| `tools/layer_lines.py RUN_DIR [...]` | how smoothly the layer lines run over a model's visible outer surface (`layer_turns_<region>` in the suite; regions in its `REGIONS`): total turning of the deformed-height contours on the surface triangles. Low = parallel bands; rings, kinks and crossing lines add to it |
| `tools/compare_gcode.py a.gcode b.gcode` | compares any two 4-axis files line by line, with tolerances |
| `tools/check_print_time.py` | checks that the pipeline's print-time estimate (`s4/print_time.py`) equals R-Theta Sim's on the fixtures in `tools/fixtures/`. Run it after changing either planner, and refresh its expected numbers from the sim (the script says how) |
| `s4_slice.py model.stl --notebook-exact` | reproduces old notebook results |
| `tools/run_reference.py map --model M.stl --out DIR --sliced planar.gcode --deformed-pkl pickle_files/deformed_M.pkl` | maps a notebook pickle with the reference port |

## How it works

1. **Deform** (`s4/fast_deform.py`, `s4/island_free.py`): tetgen, the notebook's rotation (tilt) field, then turning
   that field into a deformed shape.
2. **Slice** (`s4/cura.py`): headless CuraEngine.
3. **Map** (`s4/fast_map.py`): planar G-code mapped back to the 4-axis machine, then post-processed (axis speed
   limits, settings block, print-time estimate).
4. **Checks** (`s4/support_check.py`, `s4/quality.py`): the `[support]` and `[quality]` lines.

### The deformation (`island_free`, default)
Cura slices the deformed shape flat. Any local low point of that shape, a spot lower than everything around it that
isn't on the bed, starts printing in mid-air. `island_free` keeps the notebook's tilt field but builds the shape in
two steps:
- **fit:** first the exact least-squares fit of every tet to its target rotation (a single sparse solve), then a soft
  penalty pushes open any tet that was squashed or folded (`FIT_METHOD: penalty`). `FIT_METHOD: barrier` forbids
  folds outright instead, but it can stall on fine meshes.
- **lifting:** vertices that can't be reached from the bed by a path rising at least `ISLAND_LIFT_SLOPE` per mm get
  height targets, and the fit is solved again. Those regions then print later, growing out from where they're
  attached. How the lift holds the rest of the part still is described in the
  [CHANGELOG](../CHANGELOG.md#2026-10-01-smoother-lifted-surfaces-six-model-quality-suite). `LIFT_MODE: band`
  (opt-in; the benchy profile) fills a pit with the layer field around it instead of the flood level, so the layer
  lines cross the patch as bands instead of rings ([CHANGELOG](../CHANGELOG.md#2026-10-14-band-lift-hull-layer-bands-without-islands-lift_mode-band)).

- **pit-free field** (opt-in, `PIT_FREE_ROUNDS` > 0, between the two): the one step that changes the notebook's tilt
  field. The fit's height is, to first order, h = K^-1 b(theta) (K the mesh stiffness, b the per-tet target
  rotations' third rows); the priority flood gives each pit vertex its depth below the rim it spills over, and one
  adjoint solve gives the derivative of the total pit depth with respect to every tet's angle. That gradient, spread
  smoothly (`PIT_FREE_SMOOTH`), is the change of the field; capped steps (`PIT_FREE_MAX_STEP`), one fit per round,
  until no pit is left. Both implementations run it (it is inside `island_free`).

`DEFORMATION_METHOD: notebook` uses the notebook's own least-squares solve instead. Before/after numbers are in the
[CHANGELOG](../CHANGELOG.md#2026-09-27-island-free-deformation).

### Cura settings
Resolved from the 3mf the way the Cura GUI does it. Precedence is user > quality_changes > quality > material >
definition_changes > definition, and the extruder stack falls back to the global stack. `=expressions`, `resolve` and
`limit_to_extruder` are handled too. The result goes to CuraEngine as a generated `.def.json`, because about 700 `-s`
flags would exceed the Windows command-line limit.

### Two implementations
`--impl fast` and `--impl reference` give byte-identical output (rule 1). The speed-ups in the fast one are listed in
the [CHANGELOG](../CHANGELOG.md#2026-09-25-headless-pipeline).

### The support check
The check maps the toolpath into real space and walks it in print order. A point is **grounded** if the bed or
grounded plastic from an earlier layer is within 1 mm, or if it is within 5 mm along the same extrusion line of such
a point (normal bridging and overhang reach). Everything else is **ungrounded**: a floating island, plus every layer
stacked on top of it. The finer, local view counts points with nothing printed earlier within 1 mm, and sorts the
runs into island, cantilever and bridge by what their ends connect to. `check_planar()` runs the same check on
Cura's planar toolpath, to separate the shape's problems from the mapping's.

## Where the time goes

Every run writes `build/<model>/timings_fast.json` (the stage tree, `s4/timing.py`); the log ends with the same
tree. Seconds per stage with the same planar G-code (`--sliced-gcode`, no CuraEngine), default settings, before
(c348a0c) and after the 2026-10-14 speed-ups, run back to back on an otherwise idle 24-thread PC (the same PC under
load from other jobs is up to 3x slower):

| stage | dino (4k tets) | B axis mount 5mm (21k) | benchy upsidedown tilted (33k) | axis model w markers2 (50k) |
|---|---|---|---|---|
| mesh setup | 0.5 -> 0.4 | 2.9 -> 1.5 | 6.0 -> 1.9 | 7.8 -> 2.9 |
| rotation field (TRF least squares) | 2.6 | 6.8 | 12.1 | 12.9 |
| deformation (`island_free` fit + lift) | 0.2 | 5.1 | 10.1 | 19.6 |
| map, incl. axis speed limits | 6.8 | 8.0 | 8.8 -> 7.9 | 15.8 -> 14.7 |
| support check | 9.8 -> 1.1 | 10.2 -> 1.2 | 9.8 -> 1.0 | 23.8 -> 1.8 |
| poles + print-time estimate, after the support check | 3.2 | 2.6 -> 3.4 | 2.7 -> 3.3 | 5.6 -> 6.6 |
| whole run | 25.9 -> 16.5 | 38.8 -> 27.3 | 53.3 -> 36.8 | 90.1 -> 61.3 |

(Rows without an arrow: unchanged within the noise. The poles check and the print-time estimate now run in a worker
process during the support check; that row is the wait after it, and before it left out the untimed poles check.)
What is left is mostly the two solvers: the rotation field (scipy's TRF with LSMR; its time goes into sparse
mat-vecs and OpenBLAS's threaded dot products, not Python) and the `island_free` L-BFGS. Denser meshes
(`TET_MAX_VOLUME`) cost a little more than in proportion: the benchy at 0.3 mm^3 has 102k tets (3.1x) and its deform
took 97.5 s instead of 23.3 s (4.2x): rotation field 48.6 s (4.4x), deformation 40.7 s (4.3x), mesh setup 7.6 s
(5.2 s of it the neighbour lists, still Python sets; before the speed-ups mesh setup was 25.6 s there). The map,
the support check and the print-time estimate scale with the G-code (points), not with the mesh.

**Speed-ups must not change a byte** (rule 1 and 3). What is safe: per-element work in parallel (each tet, row or
query point computed on its own, in the same operation order), vectorising a loop with the same float operations
(a row-wise `np.linalg.norm(a, axis=1)` gives the same values as 1-D norms; numpy's `sin`/`cos` give the same values
on arrays as on scalars), caching a result that would be recomputed identically, and running independent stages at
the same time. What is not: any reduction whose summation order changes (parallel sums, BLAS dot products with a
different thread count, `np.add.at` -> `bincount`), different tie-breaking in searches (that is why the Dijkstra is a
replay of networkx's, heap order included), and Python's `x ** 2`, which is not always `x * x`. Check with the gates
and by byte-comparing the final G-code and `deformed_points.npy` of several models on the same planar G-code.

## Modules

| module | role |
|---|---|
| `s4_slice.py` | command-line entry point; parses options and calls `s4/pipeline.py` |
| `s4/pipeline.py` | runs the stages in order (deform, slice, map, checks) for either implementation; writes `build/` and the final G-code |
| `s4/profile.py` | build profiles: finds and loads the JSON, applies `--set` / `--cura-set`, rejects typos, writes `params_used.json` |
| `s4/params.py` | the default `deform` settings and the multi-iteration schedule |
| `s4/meshio_s4.py` | mesh loading, tetrahedralisation and VTK helpers shared by the fast code |
| `s4/fast_deform.py` | tilt field (notebook cells 2-11) and the notebook-method deformation, numerically equal to the reference; the cell graph's Dijkstra is a numba replay of networkx's |
| `s4/fast_lsq.py` | numba-accelerated sparse least-squares for the notebook-method solve, bit-identical to scipy's |
| `s4/island_free.py` | the `island_free` method: the fit-then-lift driver, the priority flood that finds pits, and the lift targets |
| `s4/island_free_solver.py` | the fit solver behind it: the penalty and barrier problems and their L-BFGS minimiser |
| `s4/geometry.py` | the machine geometry shared by the deformation, mapping and checks: per-tet target rotations about the tangential axis, and the nozzle-tip position from C/X/Z/B |
| `s4/cura.py` | headless CuraEngine: resolves the 3mf setting stack, writes the `.def.json`, runs the engine |
| `s4/fast_map.py` | maps Cura's planar G-code back to 4 axes (B, C, X, Z, E) and writes the final G-code |
| `s4/feed_limits.py` | post-processing: stretches each move's G93 time so no axis exceeds its speed limit |
| `s4/sim_header.py` | builds the settings block at the top of the G-code, for R-Theta Sim |
| `s4/print_time.py` | print-time estimate, a port of R-Theta Sim's parser and planner |
| `s4/support_check.py` | the `[support]` check (ungrounded, island, cantilever, bridge); reuses the deform's mesh and the mapper's mapped points, earlier-layer queries on a numba hash grid |
| `s4/quality.py` | the `[quality]` check (poles, extrusion along the nozzle axis) and the path-roughness measure |
| `s4/reference.py` | the golden reference: a line-for-line port of `main.ipynb` |
| `s4/timing.py` | nested stage timer used by both implementations |
| `tools/` | the gates and the suite (above), plus `tools/fixtures/` for the print-time check |
| `sim/r-theta-simulator.html` | R-Theta Sim, a single-file viewer ([how to use it](../sim/README.md)) |
