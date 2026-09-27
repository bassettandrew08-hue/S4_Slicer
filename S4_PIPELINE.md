# S4 headless pipeline

One command, STL in, 4-axis G-code out (C/X/Z/B/E, `G93` inverse-time feed, `M83` relative E).

## Quick start

Run these from the `S4_Slicer` folder.

End to end: deform, slice with CuraEngine (no Cura window), map back to 4 axes:
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl"
```
The output goes to `output_gcode/<model>.gcode`. Like the notebook, this **overwrites** any existing file with that
name. Intermediate files go to `build/<model>/`. A timing table for each stage is printed at the end.

Choose the output file:
```
venv\Scripts\python s4_slice.py "input_models/benchy upsidedown tilted.stl" -o my_benchy.gcode
```

Run the multi-iteration benchy recipe (the notebook's cell 4 → 7 → 9 loop):
```
venv\Scripts\python s4_slice.py "input_models/benchy upsidedown tilted.stl" --params params/benchy_upsidedown_tilted_recipe.json
```

Change a Cura setting for one run (repeat `--cura-set` for more):
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --cura-set layer_height=0.1
```

Use a different Cura profile (save it from Cura as a project file):
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --cura-config my_profile.3mf
```

Check that the fast code still matches the original notebook code exactly (run this after any code change):
```
venv\Scripts\python tools\check_equivalence.py "input_models/pi 3mm.stl"
```

Reproduce the notebook's output exactly, including its retraction behaviour (see Known issues):
```
venv\Scripts\python s4_slice.py "input_models/pi 3mm.stl" --notebook-exact
```

To see every option:
```
venv\Scripts\python s4_slice.py --help
```

## How it works

Stages (timed and printed at the end of every run, also saved to `build/<model>/timings_<impl>.json`):

1. **deform**: tetgen, rotation field, deformation solve (the notebook's cells 2–11). Parameters default to
   the `main.ipynb` values that produced the verified benchy output (`s4/params.py`: W=30, MAX_OVERHANG=30,
   ROTATION_MULTIPLIER=2, ...). Override them with `--params my.json`, either as a flat dict or as an iteration
   schedule that replays the notebook's "cell 4 → cell 7 → cell 9 → cell 4 ..." loop. Each iteration inherits the
   previous iteration's values:
   ```json
   {"iterations": [{"NEIGHBOUR_LOSS_WEIGHT": 100, "MAX_OVERHANG": 5, "ROTATION_MULTIPLIER": 1,
                    "SET_INITIAL_ROTATION_TO_ZERO": true},
                   {"NEIGHBOUR_LOSS_WEIGHT": 50}]}
   ```
   `params/benchy_upsidedown_tilted_recipe.json` is the 5-iteration recipe from notebook cell 5.
2. **slice**: headless CuraEngine (`s4/cura.py`), with settings taken from `cura_config.3mf`.
3. **map**: planar G-code mapped back to the 4-axis machine (cells 15–18).

Intermediate files go to `build/<model>/`: the deformed STL, the planar G-code, `cura_engine.log`, and `deformed_points.npy`.

## Cura settings

`s4/cura.py` resolves the 3mf's setting stack the way the Cura GUI does. Precedence is user > quality_changes > quality >
material > definition_changes > definition, with the extruder stack falling back to the global stack. `=expressions` are
evaluated, and `resolve` and `limit_to_extruder` are handled. It then enforces these settings, which the S4 mapper
depends on:

| setting | value | why |
|---|---|---|
| machine_center_is_zero | True | origin at build-plate centre |
| retraction_hop_enabled | False | no Z-hop |
| machine_gcode_flavor | RepRap (RepRap) | RepRap flavour (the printer firmware) |
| relative_extrusion | True | the mapper scales each move's E value |
| center_object | True | model centred on the plate ("Arrange All") |

Headless runs have no post-processing scripts. The resolved global settings go to CuraEngine as a generated
`.def.json` (`-j`), because about 700 `-s` flags would exceed the Windows command-line limit. Extruder values that differ
from the global ones go as `-e0 -s`. The mapper's retraction detection uses the resolved `retraction_amount`, which
is 1.0 mm in this 3mf, the same as the notebook's hard-coded `RETRACTION_LENGTH`.

Note: `cura_config.3mf` is not the profile the GUI used for the existing `input_gcode/` files. The benchy there was
sliced with layer height 0.1, **absolute** extrusion (`M82`) and 6.5 mm retraction. Absolute E breaks the mapper, which
treats every E as relative: the notebook's `output_gcode/benchy upsidedown tilted.gcode` has per-segment values like
`E1794`. To use a different profile, save it from Cura as a project (`--cura-config other.3mf`) or use `--cura-set`.

## Implementations and the equivalence gate

* `--impl reference` (`s4/reference.py`): a line-for-line port of the notebook, with only globals and plotting removed.
  It is verified **byte-identical** to the notebook's own `output_gcode/benchy upsidedown tilted.gcode`, and its
  deformed mesh is **bit-identical** to `pickle_files/deformed_benchy upsidedown tilted.pkl`.
* `--impl fast` (default) (`s4/fast_deform.py`, `s4/fast_map.py`, `s4/fast_lsq.py`): same results, faster.

```
venv\Scripts\python tools\check_equivalence.py "input_models/pi 3mm.stl"
venv\Scripts\python tools\compare_gcode.py a.gcode b.gcode      # any two 4-axis files, with tolerances
```

`check_equivalence.py` runs both implementations on the same model and the same planar G-code. It compares the
deformed vertices, the deformed STL and the final G-code, and exits non-zero on failure. Run it after any change to
`s4/`. The current fast path is not just within tolerance: its output is byte-identical.

What makes the fast path fast. Every change preserves the output:

* neighbours come from connectivity via numpy, then go through Python sets in VTK/pyvista's insertion order, so every
  neighbour list keeps its exact order (this matters for Dijkstra tie-breaking and least-squares row order)
* the unused `shared_vertices` loop is dropped; the cell↔face maps use hashing instead of an O(points × surface
  points) search
* Dijkstra runs once instead of three identical times
* the initial-rotation smoothing runs once. The notebook's 30 passes each recompute from the same unsmoothed field,
  so pass 30 equals pass 1
* least-squares Jacobians are assembled straight into CSR, with the same float32 values and pattern `lil_matrix`
  produced
* scipy TRF's sparse mat-vecs run in a row-parallel numba kernel with the same summation order (bit-identical)
* the mapper uses batched LAPACK for the Kabsch fits, volumes and barycentrics, and a regex G-code reader instead of
  pygcode. `find_closest_cell` only runs where it's needed. `calculate_tet_attributes(deformed)` is skipped, since only
  connectivity, points and cell centres are used.

## Fixed: filament "sticks" at every retraction (on by default)

The notebook marks points between a 1 mm retract (`E-1`) and the matching unretract (`E+1`) as travelling, and it
writes travelling points 1 mm further out along the tool axis (`z_hop`). So each retract came out as "lift 1 mm while
retracting", and each unretract as "plunge 1 mm while extruding 1 mm of filament", at `G94 F20000`. That's about 100×
a normal segment's E on a vertical path, which shows up as a stick of filament.

Old GUI slices used 6.5 mm retraction, which never matched the notebook's 1.0 constant, so this never fired. The 3mf's
1 mm retraction triggers it at almost every travel (about 940 times on the pi).

Now the pipeline retracts in place, then lifts, travels, lowers, and unretracts in place. The E-only lines are written
as `G94` + `F<Cura's retraction feed>` (F3600 here), then `G93`. This is option `SPLIT_RETRACTIONS` in the mapper,
and both implementations apply it identically. `--notebook-exact` turns it off and reproduces the old output byte for
byte. Note that CuraEngine itself is not deterministic: two slices of the same STL differ by about 0.01 mm on a few
hundred lines. So compare outputs using the same planar G-code (`--sliced-gcode`), as `check_equivalence.py` does.

## Fixed: blotchy extrusion multiplier (on by default)

The notebook scales every segment's E by the volume ratio of the tet it falls in. That's right on average, since it
conserves plastic when layers get stretched or squeezed, but it's constant per tet. So flow jumped wherever a path
crossed a tet boundary: on the pi, 3,366 point-to-point jumps of more than 25%. These show up as sharp red/blue
triangles in a plastic-per-mm view.

The ratio is now computed per vertex (undeformed ÷ deformed volume of the tets around it) and interpolated
barycentrically along the path, the same way position and rotation already are. It is then clamped to 0.5×–2×.
On the pi: 84 jumps over 25%, 99% of points within 0.84–1.36×, total plastic within 0.2% of before. The mapper
options are `SMOOTH_EXTRUSION_MULTIPLIER` and `EXTRUSION_MULTIPLIER_RANGE`.

## Fixed: floating blob from Cura's start code (on by default)

Cura's start code primes the nozzle (`G1 F200 E3`) at the park position (0, 0, 20). The mapper treats that point
like any other, so it became a 3 mm blob placed in mid-air inside the part's volume. The headless pipeline now sets
the start code to `G28 ; home`; the S4 header already primes at home (`G1 E10`). A `--cura-set machine_start_gcode=...`
override still wins.

`--notebook-exact` disables all three fixes (this one, the multiplier fix, and the retraction fix above).

## Support check: mid-air extrusion report

Every run ends with a `[support]` line (`s4/support_check.py`, skip with `--no-support-check`). It maps the planar
toolpath into real space, walks it in print order, and counts extruded points with nothing printed earlier (and not
the bed) within 1 mm. This is the slicer's own answer to "will this float?", independent of any simulator.

Floating runs are classified by what their ends connect to:
- **bridge**: supported extrusion at both ends. This is normal FDM bridging. Sparse gyroid infill does it even in
  a flat print: a plain 30 mm cube shows 1.3% "floating" infill with this 1 mm radius.
- **cantilever**: supported at one end only.
- **island**: connected to nothing. This is the only truly unprintable case.

Measured with default parameters:

| | pi | benchy |
|---|---|---|
| floating overall (flat print → S4) | 2.2% → 0.42% | 0.63% → 0.51% |
| islands | 3 runs, 4 mm (outer wall) | 25 runs, 21 mm |
| cantilevers | 9 runs, 31 mm | 37 runs, 106 mm |
| bridges | 49 runs, 193 mm, mostly infill | 189 runs, 488 mm |

The mapping is faithful: Cura's planar toolpath already floats in the *deformed* shape, because the deformation
doesn't make a flat 180° bridge underside (the pi's crossbar, z = 10 mm over an 11 mm gap) fully printable. So the
crossbar underside gets printed while the legs are still short of it. Tuning (MAX_OVERHANG, multiplier, weight,
iterations, part offset) moved the pi between 0.18% and 0.77%. None of it reached zero, and the best settings
drive B to its −130° limit.

## Known issues found (not changed; they would change the output)

* **The deformation solve is quartic and unconverged.** Its residual is `||N V - R N V0||²`, which least squares then
  squares again. It stops at `max_nfev=1000` (benchy: first-order optimality 1.55). The result is chaotic: nudging the
  start point by 1e-9 mm moves vertices by up to 0.6 mm. So it can only be reproduced bit-for-bit, and a GPU port can't
  match it within a small tolerance. The exact *linear* least-squares solve takes 0.2 s, and even scores lower on the
  quartic objective, but it moves benchy vertices by up to 13 mm.
* **The rotation smoothing is quartic** (`W·Δ²` residuals, so the penalty is `W²Δ⁴`). Left unchanged, as requested.
* `INITIAL_ROTATION_FIELD_SMOOTHING=30` really does one smoothing pass (see above). It was probably meant to be iterative.
* Because the deformation solve is chaotic, results can differ at the sub-mm level across machines or BLAS thread
  counts, because OpenBLAS's `ddot` summation order depends on the thread count.
