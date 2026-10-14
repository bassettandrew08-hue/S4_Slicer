# Changelog

Changes in this fork relative to [jyjblrd/S4_Slicer](https://github.com/jyjblrd/S4_Slicer). Newest first. How to use the
pipeline is in [S4_PIPELINE.md](S4_PIPELINE.md); which setting turns on which fix is in
[docs/SETTINGS.md](docs/SETTINGS.md#differences-from-the-notebook). Measurements are on the
`benchy upsidedown tilted` (33,156 tets) and `pi 3mm` models with default settings, unless noted.

**Entries**
- [2026-10-14: Faster pipeline, same output; denser meshes (TET_MAX_VOLUME)](#2026-10-14-faster-pipeline-same-output-denser-meshes-tet_max_volume)
- [2026-10-14: The benchy's choppy hull; mesh density setting](#2026-10-14-the-benchys-choppy-hull-mesh-density-setting)
- [2026-10-14: Where the path roughness comes from; opt-in smoother mapping](#2026-10-14-where-the-path-roughness-comes-from-opt-in-smoother-mapping)
- [2026-10-07: Meshes with CAD-export junk, support check on large prints](#2026-10-07-meshes-with-cad-export-junk-support-check-on-large-prints)
- [2026-10-01: Smoother lifted surfaces, six-model quality suite](#2026-10-01-smoother-lifted-surfaces-six-model-quality-suite)
- [2026-10-01: Bug fixes, cleanup, path-roughness view](#2026-10-01-bug-fixes-cleanup-path-roughness-view)
- [2026-09-28: R-Theta Sim in the repo, print times in sync](#2026-09-28-r-theta-sim-in-the-repo-print-times-in-sync)
- [2026-09-28: Realistic axis speeds](#2026-09-28-realistic-axis-speeds)
- [2026-09-28: Nozzle tilt on every model (Squirtle had none)](#2026-09-28-nozzle-tilt-on-every-model-squirtle-had-none)
- [2026-09-27: Settings block for R-Theta Sim](#2026-09-27-settings-block-for-r-theta-sim)
- [2026-09-27: No more "poles" at the base](#2026-09-27-no-more-poles-at-the-base)
- [2026-09-27: Island-free deformation](#2026-09-27-island-free-deformation)
- [2026-09-26: Output fixes, build profiles, docs](#2026-09-26-output-fixes-build-profiles-docs)
- [2026-09-25: Headless pipeline](#2026-09-25-headless-pipeline)

---

## 2026-10-14: Faster pipeline, same output; denser meshes (TET_MAX_VOLUME)

Output unchanged: the final G-code and `deformed_points.npy` of benchy upsidedown tilted, B axis mount 5mm, dino and
axis model w markers2 are byte-identical to before (same planar G-code; all but the `pipeline_version` /
`planar_gcode` header lines), and the gates pass (notebook-exact, equivalence twice, print-time parity).

- **Mesh setup:** the cell graph's edge lengths in one row-wise norm (the per-edge loop of 1-D norms took 13 s on
  the benchy); networkx's graph and `multi_source_dijkstra` replaced by a numba replay of the same search (heap
  order by distance then push count, networkx's neighbour order, the same additions), so distances, closest bottom
  cells and the in-air test are identical. Benchy 6.0 -> 1.9 s, axis model 7.8 -> 2.9 s.
- **Support check** (9.8 -> 1.0 s on the benchy, 23.8 -> 1.8 s on the axis model): reuses the deform's tet mesh (no
  second TetGen run) and the mapper's parse and mapping of the planar points (the same numbers); layers and types
  from the parse's line numbers; the "anything printed earlier within 1 mm" queries on one static numba hash grid
  instead of KD-trees rebuilt every layer (same distance test, same answers on all four models); the distance to the
  surface only for points outside the part.
- **Map:** the planar G-code parse vectorised (same float operations, ~3x faster) and cached for the support
  check; axis speed limits compute all move times at once.
- **Poles check and print-time estimate** run in a worker process while the support check runs (`S4_NO_WORKER=1`
  runs them in-process). The end of the run is timed now too (`poles + print-time estimate`, `settings header`).
- **`TET_MAX_VOLUME`** (deform, default null = TetGen's default mesh, unchanged): largest tet volume in mm^3, for
  denser meshes. Benchy: 0.3 -> 102k tets (3.1x), deform 97.5 s instead of 23.3 s (rotation field 48.6 s,
  deformation 40.7 s: a little more than in proportion). The header says `tetrahedralize = fixedvolume, maxvolume ...` when it is set.

Whole runs (same planar G-code, back to back on an idle PC): dino 25.9 -> 16.5 s, B axis mount 38.8 -> 27.3 s,
benchy 53.3 -> 36.8 s, axis model 90.1 -> 61.3 s. Per stage: [DEVELOPING.md](docs/DEVELOPING.md#where-the-time-goes). Not faster yet: the rotation-field
least squares and the `island_free` solve, now most of the time. A numba port of scipy's LSMR (bit-identical, it
calls numpy's own BLAS dot) was tried and dropped: the time is in OpenBLAS's threaded dot products (~0.5 ms each
here) and the sparse mat-vecs, not in Python, so it saved nothing.

---

## 2026-10-14: The benchy's choppy hull; mesh density setting

### Diagnosis
The chop on the benchy's hull (sides of the hull, above the cabin's front and back ends, z ~28-42 in the part's
frame, both sides) is where the island lift reshaped the surface. Every layer line on the outer surface is a contour
of the deformed height; drawn on the surface's triangles (scratch render, and in R-Theta Sim), the hull's lines run
in smooth parallel bands everywhere except in the lifted patch, where they turn into wide, wandering, polygonal
loops that cross the neighbouring bands. Two parts to it:
- The lift fills a pit to the "flood" level, a surface rising at `ISLAND_LIFT_SLOPE` from the pit's rim, so the hull
  there becomes nearly parallel to the layers (grazing): lines far apart, closed rings around the high points. 921
  of the 953 pit vertices are surface vertices: the pit is the hull's skin.
- Where the layers graze the surface, mesh-scale height noise (the flood is a graph-path cone; the lift rounds add
  targets round by round) moves the lines a lot: they follow the triangle edges.
The fit alone (no lift) draws the hull smoothly: hull `layer_turns` 384 vs 538 with the lift, 200 for a flat print.

### Added
- **`tools/layer_lines.py`** and the suite metric **`layer_turns_hull`** (benchy only for now): the total turning
  of the layer lines over the visible hull exterior (z > 26.5), in full turns, computed from the deformed mesh on
  the surface triangles. A closed ring counts 1; lines that kink at tet faces or wander where layers graze count a
  lot. It follows what the sim shows, and is about mesh-independent (a flat print: 200 on every mesh tried).
  `layer_kink_turns_hull` counts only kinks over 15 deg at one edge.
- **`TET_MAX_VOLUME`** (deform, default 0 = off): caps every tet's volume (TetGen `-a`), a finer mesh. Used by the
  deformation, the reference port and the support check alike, so the support check re-tetrahedralizes the same mesh.

### Changed
- **`params/benchy upsidedown tilted.json`** (new; benchy only): `MAX_OVERHANG` 20, `ISLAND_LIFT_SLOPE` 0.2,
  `LIFT_HOLD_FALLOFF` 5. More tilt leaves a smaller pit for the lift (570 lifted vertices instead of 977), and the
  gentler lift with a soft-edged hold puts fewer crossing lines into the hull. The sim shows clean parallel bands
  over most of the hull, with concentric rings (smooth, not crossing) where the lift remains.

| benchy | default (branch) | profile |
|---|---|---|
| hull `layer_turns` (kinks) | 538 (301) | 411 |
| ungrounded total (walls/skin) | 41.4 (21.9) mm | 19.8 (1.8) mm |
| island walls / cantilever walls | 10.5 / 1.0 mm | 2.6 / 6.9 mm |
| zigzag / surface_rough_pct / waviness | 4629 / 1.06 / 674 | 3937 / 0.92 / 185 |
| tilt | 49.4 deg | 56.5 deg |
| print time | 151.3 min | 160.7 min |

The one thing that gets worse is cantilevered wall (1.0 -> 6.9 mm); unsupported and island walls drop much more.

### Tried and not adopted (hull `layer_turns`, deform only, same rotation field; default 538)
| change | layer_turns | why not |
|---|---|---|
| `TET_MAX_VOLUME` 1 / 0.5 / 0.25 (45k / 67k / 122k tets) | 524 / 530 / 608 | no better; deform 94 / ~500 / 308 s vs 91 s (lift solves dominate) |
| `ISLAND_LIFT_SLOPE` 0.3 / 0.2 | 470 / 430 | in the profile (0.2) |
| `LIFT_HOLD_FALLOFF` 5 / 10 | 494 / 475 | in the profile (5) |
| `LIFT_ITERATIONS` 300, `LIFT_HOLD` 0, two-sided lift | 538 / 532 / 526 | no effect |
| old (unanchored) lift | 624 | worse |
| smoothing the lift's displacement afterwards (20-200 passes) | 497-395 | undoes the lift: island mass 0 -> 25-100 |
| smooth upper bound of the lift targets (dilate + Gaussian, 1-4 mm) | 572-754 | lifts more area: worse |
| refit: the lift's extra rotations (smoothed over 3 mm) added to the target rotations, fit and lift again, x2 | 465 | with the profile: 413 vs 411, cantilever walls 15.9 mm |
| `MAX_OVERHANG` 20 alone / 25 + slope 0.2 | 480 / 412 | cantilever walls 8.2 / 10.5 mm |

The code for the three lift experiments is not kept.

---

## 2026-10-14: Where the path roughness comes from; opt-in smoother mapping

### Diagnosis (`tools/roughness_breakdown.py`, six models, default settings)
Every rough bead segment (the suite's `zigzag`) traced back to the planar point it came from:

| model | zigzag | rough in Cura's planar path too | planar corner + mapping | tet face / rotation split / other mapping | flat print of the undeformed model, same spots rough |
|---|---|---|---|---|---|
| pi | 213 | 195 | 18 | 0 | 14 |
| benchy | 5228 | 4163 | 800 | 265 | 2957 |
| Squirtle | 2134 | 1841 | 231 | 62 | 1341 |
| dino | 3163 | 3019 | 116 | 28 | 2074 |
| z mount | 1387 | 1183 | 157 | 47 | 336 |
| B axis mount | 3542 | 3027 | 454 | 61 | 983 |
| **all** | **15667** | **13428 (86%)** | **1776 (11%)** | **463 (3%)** | **7705 (49%)** |

- About half of the count is the model's own geometry as Cura draws it: small fillets drawn with 2-3 segments
  turning 30-45 deg each, skin outlines, infill ends. A flat print of the undeformed STL is rough at the same spots
  (flat-print totals: benchy 4679, Squirtle 2310, dino 2763, z mount 1537, B axis 1244, pi 10).
- Most of the rest is also already in Cura's planar path, but only because of the deformed shape (B axis mount and
  pi above all): contours where a layer grazes a nearly flat deformed surface, tiny closed loops at the tops of bumps.
- The mapping itself adds little above 30 deg. Below it, it makes the walls wavy: where the planar path runs
  straight, 5% of the joints of the final path turn more than 5 deg (z mount: 9837 of 191240). Two causes:
  the per-tet linear map bends lines at tet faces, and the machine moves C and X linearly, so near the turntable
  axis a move that sweeps several degrees of C is drawn as a curve (almost all turns > 10 deg are within 5 mm of
  the axis).
- Cura's `meshfix_maximum_deviation` 0.01 (profile: 0.025) draws curves with more, gentler segments: zigzag
  benchy 5228 -> 4853, Squirtle 2134 -> 1659, dino 3163 -> 2539, z mount 1387 -> 1185, B axis 3542 -> 2442, but
  pi 213 -> 267 (bridge underside 31 -> 41) and z mount ungrounded walls 0.3 -> 2.7 mm. Not adopted; it is a
  `--cura-set` away for the models it helps. 0.005 and 0.05 are both worse than 0.01.

### Added (opt-in, off by default: no output change)
- **`MAP_SMOOTH`** (mm) and **`MAP_SMOOTH_MAX`**: the displacement at every printed point is re-fitted by a local
  linear fit over its neighbours on the same printed run (`s4/fast_map.py smooth_mapping`). Affine maps, corners
  and curves are reproduced exactly; only the bends at tet faces are blended. Points move 0.005-0.02 mm on average.
- **`MAX_C_STEP`** (deg): print moves are split so none sweeps more than this much C.
- Both are in `s4/reference.py` too (it calls the same smoothing function), so fast and reference stay
  byte-identical with them on: `tools/check_equivalence.py "input_models/pi 3mm.stl" --set MAP_SMOOTH=0.5
  --set MAX_C_STEP=3` passes. `check_equivalence.py --set KEY=VALUE` is new and goes to both runs.
- **`tools/roughness_breakdown.py`**: the breakdown above. **Suite metric `waviness`**: joints of the final path
  turning > 10 deg where Cura's planar path runs straight (what the mapping adds; uses the new `map_src.npy` in the
  work dir, written by `--impl fast`). **`tools/model_suite.py --sliced-from DIR`** reuses an earlier run's planar
  G-code, for clean comparisons of map-only settings.

### Fixed
- **`outside_mm` false positives** (`tools/shape_metrics.py`): the sign came from open3d's single-ray test, which
  flips on meshes that aren't watertight (Squirtle's STL is open and self-intersecting). A segment 1.5 mm *inside*
  Squirtle, 0.001 mm from a segment of the default run, counted as 1.5 mm outside. The sign is now a 5-ray vote;
  `surface_rough_pct` (unsigned) is unchanged.

### Results: six models, re-sliced (no `--sliced-from`), default vs `MAP_SMOOTH` 0.5 + `MAX_C_STEP` 3

| model | zigzag | surface_rough_pct | waviness | minutes |
|---|---|---|---|---|
| pi (bridge underside) | 213 (31) -> 195 (22) | 0.17 -> 0.16 | 1986 -> 191 | 40.6 -> 41.2 |
| benchy | 5228 -> 4628 | 1.22 -> 1.06 | 6539 -> 672 | 149.9 -> 151.4 |
| Squirtle | 2134 -> 2030 | 0.84 -> 0.80 | 1756 -> 111 | 82.1 -> 82.8 |
| dino | 3161 -> 3164 | 1.24 -> 1.24 | 3310 -> 78 | 98.2 -> 99.0 |
| z mount | 1389 -> 1292 | 0.42 -> 0.39 | 2516 -> 137 | 93.5 -> 94.4 |
| B axis mount | 3544 -> 3554 | 1.27 -> 1.25 | 3257 -> 231 | 117.4 -> 118.1 |

Support numbers (ungrounded, island, cantilever, bridge), poles, along-axis, tilt and `outside_mm` (0) are
identical on all six. Waviness drops 90-98%. On z mount alone (same planar G-code), joints turning > 5 / > 10 deg:
9837 / 2514 by default, 6178 / 1450 with `MAP_SMOOTH` 1.0 alone, 5771 / 1149 with `MAX_C_STEP` 3 alone, 527 / 40
with both. `MAP_SMOOTH` 1.0 instead of 0.5 is a little worse on B axis mount (zigzag 3616, surface 1.28).
A default candidate; not switched on here.

---|---|---|---|
| pi (bridge underside) | 213 (31) / 0.17 | 195 (22) / 0.16 | 200 (22) / 0.16 |
| benchy | 5228 / 1.22 | 4628 / 1.06 | 4680 / 1.08 |
| Squirtle | 2134 / 0.84 | 2030 / 0.80 | 2040 / 0.80 |
| dino | 3163 / 1.24 | 3159 / 1.24 | 3111 / 1.22 |
| z mount | 1387 / 0.42 | 1290 / 0.39 | 1306 / 0.39 |
| B axis mount | 3542 / 1.26 | 3555 / 1.25 | 3616 / 1.28 |

Wobble (joints turning > 5 / > 10 deg where the planar path is straight), z mount: 9837 / 2514 by default,
6178 / 1450 with `MAP_SMOOTH` 1.0 alone, 5771 / 1149 with `MAX_C_STEP` 3 alone, 527 / 40 with both. B axis mount
with `MAP_SMOOTH` 0.5 + `MAX_C_STEP` 3: 1757 / 231. Support numbers, poles and along-axis are unchanged (the support
check maps the planar path itself); print time +0.5-0.8 min. One problem: with `MAP_SMOOTH` 0.5 Squirtle has one
0.5 mm segment 1.5 mm outside the part (`outside_mm` 0 -> 0.5), not seen with 1.0; not investigated yet, so neither
is the default.

---

## 2026-10-07: Meshes with CAD-export junk, support check on large prints

### Fixed
- **STLs that TetGen rejects** ("The input surface mesh contain self-intersections"), e.g. `axis model w markers2`:
  the part itself was watertight, but the CAD export left zero-area sliver triangles around the markers. When
  TetGen rejects a surface, `meshio_s4.tetrahedralize` now retries once on a cleaned copy (vertices closer than
  0.001 mm merged, sliver and duplicate triangles dropped, loose flat pieces with no volume dropped) and prints
  `[mesh] ... Retrying`. Meshes that TetGen accepts as they are never reach the cleanup, so every other model's
  output is unchanged; the reference implementation uses the same call.
- **Support check out of memory or running for hours on large prints** (`axis model w markers2`: MemoryError,
  20+ GB, 28 minutes): it listed every pair of path points within 1 mm, about 2.7 billion pairs there. It now sweeps
  the layers in print order and only asks "is anything from an earlier layer within 1 mm" (a count, from KD-trees
  of the earlier points), never listing neighbours: 15 s and 0.3 GB on that model. Byte-identical reports on all six
  suite models, and a little faster. When a bad mapping piles thousands of points into one spot, the ungrounded
  regions in the report are grouped by 1 mm grid cells instead of exact pairs, and the report says so.
- **Support check: a point with more than 127 grounded neighbours could count as floating.** The old neighbour
  count was an 8-bit integer and wrapped negative; the new check doesn't count in 8 bits. The six suite models'
  numbers don't change.
- TetGen's `_skipped.face` / `_skipped.node` dumps of a rejected surface are deleted after the cleanup retry.

### Added
- **`SURFACE_SIMPLIFY_ERROR`** (deform, default 0 = off) and `params/axis model w markers2.json`, which sets it to
  1e-5. That model's marker holes are drawn with 0.04 mm segments; TetGen filled them with microscopic tets
  (399,265 tets for an 8,500-triangle surface) and the deformation tore a cluster of them apart at the hole at
  z = 11, so about 12 m of wall toolpath mapped into that 2 mm hole. The setting coarsens needlessly fine
  triangles before TetGen (quadric edge collapse with an error bound) and logs how far the surface moved:

  | axis model w markers2 | default | `SURFACE_SIMPLIFY_ERROR` 1e-5 |
  |---|---|---|
  | surface triangles / tets | 8,500 / 399,265 | 3,232 / 50,292 |
  | surface moved | - | at most 0.012 mm |
  | inverted tets | 44,679 | 818 |
  | ungrounded (walls/skin) | 12,767 mm (11,780) | 47 mm (47) |
  | island / cantilever / bridge | 27.6 / 47.3 / 215.3 mm | 23.9 / 39.5 / 103.8 mm |
  | extrusion along the nozzle axis | 163 moves | none |
  | full run | 12 min | 1.5-4 min |

  Off by default, so every other model's mesh and output are unchanged (`--notebook-exact` body identical,
  fast = reference, six-model suite unchanged).

---

## 2026-10-01: Smoother lifted surfaces, six-model quality suite

### Changed
- **New default island lift** (`s4/island_free.py`). The pi's bridge underside printed as a rough zig-zag patch.
  The cause was the lift, not the turntable axis the bridge sits on:
  - nothing held the part down while lifting, so the lift mostly moved the whole part up (pi: ~4.9 mm everywhere,
    ~1 mm of real lift under the bridge)
  - targets kept from earlier rounds, being two-sided, then pulled vertices back down: creases
  - a pit's rim got dragged up round after round

  Now the grounded vertices are anchored (`LIFT_ANCHOR` 1), only the pits are pushed up (`LIFT_ONE_SIDED`), and
  every vertex that isn't in a pit is held where the fit put it (`LIFT_HOLD` 0.2). The lift solves also get their
  own solver stiffness floor (`LIFT_PRECOND_FLOOR` 0.3); without it, Squirtle's dense cluster of tiny tets was
  torn apart and printed about 5 mm of wall as a blob. The old lift is `LIFT_ANCHOR=0 LIFT_ONE_SIDED=false
  LIFT_HOLD=0 LIFT_PRECOND_FLOOR=-1`. The part stays where it was, and the rotation field (Byrd's math) is
  unchanged.

### Results (six models, default settings; unsupported walls/skin = the part that matters structurally)

| model | unsupported, total (walls/skin) | islands | cantilevers, walls/skin | roughness | pi bridge underside |
|---|---|---|---|---|---|
| pi | 0 → 0 (0 → 0) | 0 → 0 | 0 → 0 | 290 → 213 | 64 → 31 |
| benchy | 57.3 → 41.4 (23.1 → 21.9) | 26.8 → 10.5 | 4.1 → 1.0 | 5576 → 5231 | |
| Squirtle | 4.8 → 0.6 (4.8 → 0.7) | 0 → 1.3 | 8.4 → 4.6 | 2659 → 2134 | |
| dino | 11.1 → 11.1 (5.0 → 4.9) | 8.6 → 8.4 | 25.1 → 24.5 | 3174 → 3163 | |
| z mount | 2.7 → 0.3 (3.6 → 0.3) | 4.1 → 0.8 | 2.2 → 0 | 1564 → 1389 | |
| B-axis mount | 2.7 → 1.5 (2.2 → 1.3) | 3.4 → 1.8 | 12.5 → 13.0 | 3695 → 3541 | |

Old lift vs new lift on the same code and metrics (`tools/model_suite.py`, old lift via its settings). Squeezed
volume (tets below half their volume) is about the same: benchy 3.3% → 2.3%, others within ±0.3%, none inverted.

No poles, no extrusion along the nozzle axis, tilt unchanged or slightly up. The one new item is Squirtle's
1.3 mm skin dot (two points), which is in the fit's own shape: it prints with lifting off too.

### Added
- **`tools/model_suite.py`:** slices all six models and flags regressions against a baseline run.
- **`s4/quality.py`:** the sim's path-roughness measure in Python.
- **Support check numbers:**
  - true unsupported path length (it counted 0.3 mm per point, but points are ~0.5 mm apart)
  - plastic piled into one spot (a blob) counts too
  - walls/skin vs sparse infill for every class: `ungrounded_wall_mm`, `island_wall_mm`, `cantilever_wall_mm`

### Fixed
- **`tools/check_equivalence.py`** now requires the final G-code to be byte-identical, apart from the
  `planar_gcode` line (the `--sliced-gcode` run names its planar file). Before, it only printed the result.

### Tried and not adopted (measured on all six models)
- Fading the tilt near the turntable axis: no help for the pi, and it cost Squirtle and benchy support.
- Smaller `SEG_SIZE`, or splitting paths exactly at tet faces: both made the pi bridge worse (the exact path is
  itself rough where tets are squashed; the 0.6 mm sampling hides some of it).
- A minimum tet weight in the fit: didn't close Squirtle's tear cleanly.
- Moving the part off the axis: helped the bridge, but moved the part and left cantilevers at the pi pillars.

---

## 2026-10-01: Bug fixes, cleanup, path-roughness view

### Added
- **Sim: "Path roughness" colour mode.** Each bead is coloured by how sharply the path turns at it. A turn only counts
  when a neighbouring segment also turns sharply, so a part's real corners stay neutral and zig-zag patches light up
  (e.g. a rough bridge underside).

### Fixed
- **The first retract of every print swung C up to 179° with no speed limit.** It was written as a G94 move at the
  first point's position, which the speed limiter skipped. The limiter now turns every G94 line that moves into a
  timed G93 move, using its feed and the axis limits.
- **`--sliced-gcode` assumed a 1 mm retraction.** The mapper only recognises retract/unretract moves of exactly the
  retraction length. It's now read from the planar file (its most common E-only retract).
- **Pole and along-axis checks** used a hard-coded 1 mm unretract. They now use the run's retraction length.
- **`--sliced-gcode` runs** now write `; s4: planar_gcode = ...` into the settings block, as documented.
- **A `null` speed or acceleration in a profile** crashed the run at the end. It now skips the estimate with a
  message.
- **`moves_slowed_pct`** no longer counts moves converted from G94 as slowed.
- **`tools/compare_gcode.py`** silently ignored `--tol-f-rel` values below 1e-4.
- **Sim:** opening a file with no moves left the previous print on screen while the sim already pointed at the new
  file. The tip-speed readout also measured across pauses and file changes.

### Changed
- **Code cleanup, no change in output:** unused imports, arguments and constants removed, inline imports hoisted,
  stale docstrings fixed.
- **Docs split by reader:** `S4_PIPELINE.md` is the user guide; the settings reference moved to
  `docs/SETTINGS.md`, the developer notes to `docs/DEVELOPING.md`, and R-Theta Sim's usage to `sim/README.md`.
  `README.md` gained a quick start and a map of the repository.

---

## 2026-09-28: R-Theta Sim in the repo, print times in sync

### Added
- **R-Theta Sim** is now in the repo: `sim/r-theta-simulator.html` (single file; open it in a browser).
- **The sim's print time and the pipeline's estimate now agree.** Before, the sim showed 2.2× to 3.4× the
  pipeline's number, and the ratio changed from model to model. Two things differed:
  - the machine limits: the sim's placeholder Z speed was 12 mm/s against the pipeline's 50, which hit Z-heavy
    models like the benchy hardest
  - the time model: the sim plans acceleration and corners; the pipeline added up per-move times

  The pipeline now estimates with a line-by-line port of the sim's parser and planner (`s4/print_time.py`), and
  both use the same machine defaults. `tools/check_print_time.py` checks that they give the same total: they match
  exactly on its fixtures.
- **Machine settings for the estimate:** `MAX_SPEED_E`, `MAX_ACCEL_X/Z/B/C/E`, `CORNER_SPEED`, `HOME_X/Z/B`,
  `HOME_SPEED`. They're placeholders equal to the sim's defaults, and they're written to the settings block.
- **Sim, File panel:**
  - the slicer's estimate next to the sim's time
  - a table of machine values that differ between the file and the sim, with a "Use the file's values" button
    (never applied on its own)
- **Sim settings are remembered** in the browser. Only values changed from the defaults are stored, so later
  default fixes still apply.

### Fixed (sim)
- **Every print was drawn mirrored.** The pipeline writes C = atan2(y, x). With the sim's C direction set to
  "Normal", a counter-clockwise path was drawn clockwise seen from above. The default is now "Reversed", which
  matches the S4 convention. Settings files saved before this change bring the old value back.
- **False "Past axis limit" flags:** the default B range was −100° to 40°, and the pipeline tilts to −130°. The
  benchy had 14,030 legal moves flagged. The default is now −130° to 30°, the pipeline's range.
- **Default speeds** now equal the pipeline's `MAX_SPEED_*` (X 150, Z 50 mm/s, B 180, C 360°/s).
- **Playback at 1× was slower than real time** below 10 fps, and the tip-speed readout was too low with it.
- **Messages:** removed stale references to the S4 notebook and to an `s4_fix_strays.py` that doesn't exist.
  Stretched G93 moves are now blamed on "acceleration or axis limits".
- **Offline:** without internet, the page now says three.js couldn't load instead of staying blank.

### Fixed (pipeline)
- **Newlines in settings-block values** (a profile description, a Cura override) are escaped. A raw newline would
  have ended the comment, and the rest of the line would have become G-code.

### Results

Pipeline estimate vs the sim (default settings), pi: 41:12 vs 41:09. The estimate is rounded to 0.1 min.

---

## 2026-09-28: Realistic axis speeds

### Fixed
- **Moves far too fast for the machine** (seen in R-Theta Sim). Each G93 move got the time of its planar move
  (planar length / Cura's feed), which ignores what the axes do to get the nozzle tip there:
  - a 1° tilt swings X/Z ~0.7 mm through the nozzle offset
  - a small move near the centre can need a big C rotation
  - a travel lifted over the part is much longer in real space than in the planar file

  The benchy asked for up to 21,000°/s on C while printing and 900,000°/s on travels. Now each move takes at least as
  long as every axis needs at its speed limit (`s4/feed_limits.py`, after the G-code is written, so fast and
  reference stay identical). The G94 `F20000` moves become G93 too, so every move has a definite time. Moves are
  only slowed, never sped up. `--notebook-exact` leaves the feeds alone.

### Added
- **Settings `LIMIT_AXIS_SPEEDS`, `MAX_SPEED_C` (360°/s), `MAX_SPEED_B` (180°/s), `MAX_SPEED_X` (150 mm/s),
  `MAX_SPEED_Z` (50 mm/s).** These are placeholders; set your machine's real limits. They're written to the sim
  header, along with `estimated_print_time_min` and `moves_slowed_pct`.

### Results (default settings)

| model | estimated print time | moves slowed |
|---|---|---|
| Squirtle | 44 min | 50% |
| benchy | 91 min (39 min at the old, impossible speeds) | 78% |
| pi | 24 min | 66% |
| dino | 56 min | 51% |
| z mount | 54 min | 44% |
| B-axis mount | 65 min | 47% |

Only the timing changes. The toolpaths, supports and tilts are the same as before (pi 0 mm ungrounded, benchy 57 mm).

---

## 2026-09-28: Nozzle tilt on every model (Squirtle had none)

### Fixed
- **No nozzle tilt on Squirtle.** Its B stayed within ±1.5°, so it printed essentially flat, with islands. The
  fold-free solver stalled after moving vertices about 0.3 mm. There were two causes:
  - **A bug in the solver's preconditioner** (the stiffness matrix that steers each step). It was built with each
    tet's basis-function gradients transposed: the columns of `Dm⁻¹` instead of its rows. The energy and gradient
    were right, so results were valid, just badly steered. That was harmless on simple meshes and crippling on
    Squirtle's fine, irregular one: 115k tets, 11% of them micro-tets.
  - **The hard fold barrier.** On meshes with slivers and micro-tets, a few tets pin every step.
- **New default `FIT_METHOD: "penalty"`.** Start from the exact least-squares fit (one sparse solve, which reaches the
  target tilt), then a soft penalty pushes open any tet squashed below 0.5× volume or folded. There's no step cap,
  so nothing can jam. The old method is `FIT_METHOD: "barrier"`, now with a correct preconditioner and optional
  sliver/micro-tet exemptions.
- **Defaults `FOLD_PENALTY_DET: 0.5` (was 0.2) and `LIFT_WEIGHT: 50` (was 5).** At 0.2, thin features squashed flat
  and the pi extruded along the nozzle axis (278 segments). The lift weight goes up with the penalty so lifting still
  wins (at 5, Squirtle kept an unsupported tower).
- **Poles on the dino** (2, one of them 21 mm long). A travel ended just outside the part and mapped 0.07 mm below the
  bed. `write_gcode` drops points below z = 0, so the travel down to the re-entry point disappeared (and the unretract
  with it), and the next print move extruded straight down from the travel height. With `SAFE_TRAVEL_TRANSITIONS`
  (default; off with `--notebook-exact`), mapped points below the bed are now clamped to z = 0.
- **`[quality]` pole line numbers** now count the settings block, so they match the final file.
- **The bed face always counts as supported.** Island lifting only treated the lowest deformed points as standing on
  the bed. The exact fit tilts the base with the region above it, so most of a part's bed face ended up higher and
  was lifted as if it floated. Now vertices on the real bed are grounded wherever they land (`BED_TOL`, 0.3 mm).
  `BED_PIN_WEIGHT` (off by default) can also hold the bed face flat, but that fights the tilt near the base: the
  benchy got 1475 mm of unsupported extrusion, the B-axis mount lost 5° of tilt.
- **`[quality]` no longer flags sideways printing.** It used to count every near-vertical extruding move as "steep".
  With the nozzle tilted toward −90°, a vertical move is just printing along a curved layer, which is the point of
  S4. It now flags extrusion along the nozzle's own axis (within ~27°), which pushes into or pulls out of the bead.

### Added
- **`[deform]` now reports the tilt achieved against the tilt aimed for:** volume-weighted rotation about the B axis,
  where more than 20° is wanted. It warns if the tilt falls short, so a stalled deformation can't pass silently again.
  It also reports the folded volume share.

### Results (default settings)

| model | tilt achieved / aimed | ungrounded extrusion | extruding along the nozzle axis |
|---|---|---|---|
| Squirtle | 24.8° / 27.7° (B was within ±1.5°) | 5 mm (was 62 mm) | 0 |
| benchy | 48.7° / 52.0° | 57 mm | 0 |
| pi | 37.5° / 39.8° | 0 mm | 0 (278 at the old defaults) |
| dino | 22.3° / 30.1° | 11 mm | 0 |
| z mount | 27.0° / 32.2° | 3 mm | 0 |
| B-axis mount | 28.9° / 33.4° | 3 mm | 0 |

Large tilts are expected: parts like the dino's nose print nearly sideways (B toward −90°), in rings growing outward.
No model has poles. Squirtle's deformation takes about 85 s.

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
- `;SETTING_3` is written exactly as Cura writes it: the JSON is escaped (each `\` becomes `\\`, each newline
  becomes `\n`) before it's split into 80-character comment lines. A reader has to undo that escaping before parsing the JSON, as Cura's
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
- Where the islands come from: Cura slices the deformed shape flat, so every local low point of it (not on the bed)
  starts printing in mid-air. These low points match the toolpath islands exactly; the benchy's two islands are its
  two largest.
- Why the notebook leaves them: its solve stops at 1000 evaluations, unconverged. It misses its own target rotations
  by 9° median (26° at the 95th percentile) and folds 202 tets. Near the islands the target tilt is also too small
  (flat downward faces need about 60°; the smoothing settles at about 43°).
- Joshua's original has them too: his committed pi toolpath has 17 mm of ungrounded extrusion, about the same as the
  notebook method here (measured on his planar toolpath, since tetgen versions differ).

### Tried and not adopted
- Parameter tuning: the pi moved between 0.18% and 0.77% floating, never to zero, and B hit its −130° limit.
- Moving the part off-centre: no help.
- An exact linear deformation solve: better fit, but more floating extrusion (2,000 mm).
- Quadratic rotation smoothing, per-cell overhang feedback: worse.
- Lifting heights directly: inverted thousands of tets.
- Pinning the bottom face flat: jammed the solver.

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

### Fixed (on by default; `--notebook-exact` restores the notebook's behaviour, [which setting is which](docs/SETTINGS.md#differences-from-the-notebook))
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
- The floating plastic over the pi's bridge comes from the slicer, not the simulator: Cura's own planar toolpath
  already floats in the deformed shape.
- Sparse gyroid infill "floats" by nature: a plain 30 mm cube shows 1.3% under a 1 mm-radius check. Bridges anchored
  at both ends are normal.
- CuraEngine isn't deterministic: two slices of the same STL differ by about 0.01 mm on a few hundred lines, even
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
- The notebook's deformation solve is quartic, unconverged and chaotic: nudging the start point by 1e-9 mm moves
  vertices by up to 0.6 mm. It can only be reproduced bit-for-bit (a GPU port couldn't match it), so it was sped up
  on the CPU instead.
- The rotation smoothing is also quartic (`W·Δ²` residuals). Left unchanged.
- `INITIAL_ROTATION_FIELD_SMOOTHING = 30` does one pass: each of the 30 recomputes from the same unsmoothed field.
- The old GUI-sliced benchy used absolute extrusion (M82). The mapper treats E as relative, so that output had
  per-segment values like `E1794`.
