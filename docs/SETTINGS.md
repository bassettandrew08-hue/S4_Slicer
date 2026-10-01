# Settings reference

Every setting the pipeline understands, with its default. For how to *use* settings (`--set`, profile files,
variants), see the tutorial in [S4_PIPELINE.md](../S4_PIPELINE.md#3-tutorial-changing-settings). Measurements and the
reasons behind the defaults are in the [CHANGELOG](../CHANGELOG.md), not here.

**Contents**
- [Profile file format](#profile-file-format)
- [`deform`: how the part is warped](#deform-how-the-part-is-warped)
- [`map`: the machine and the final 4-axis moves](#map-the-machine-and-the-final-4-axis-moves)
- [`cura`: slicing](#cura-slicing)
- [Differences from the notebook](#differences-from-the-notebook) (which setting turns on which fix)

## Profile file format

Every section and key is optional; anything left out keeps its default.
```json
{
  "description": "free text",
  "deform": {"MAX_OVERHANG": 10, "PART_OFFSET": [0, 0, 0],
             "iterations": [{"NEIGHBOUR_LOSS_WEIGHT": 100}, {"NEIGHBOUR_LOSS_WEIGHT": 50}]},
  "map":    {"NOZZLE_OFFSET": 42, "MIN_ROTATION": -120},
  "cura":   {"config": "cura_config.3mf", "set": {"layer_height": 0.1}, "strip_start_prime": true}
}
```
- **Iterations:** `iterations` replays the notebook's "run cell 4 → 7 → 9 again" loop. Each entry lists only what
  changes, and inherits everything else from the previous iteration.
- **Degrees:** the rotation-limit keys the notebook kept in radians can be written in degrees with a `_DEG` suffix:
  `MAX_POS_ROTATION_DEG`, `MAX_NEG_ROTATION_DEG`, `ROTATION_MAX_DELTA_DEG`.
- **Old files:** a flat dict, or `{"iterations": [...]}`, still loads and is read as the `deform` section.

## `deform`: how the part is warped

Notebook cells 2-7. The first table is the notebook's rotation (tilt) field; the second builds the deformed shape
from it.

### Tilt field
| setting | default | what it does |
|---|---|---|
| `MAX_OVERHANG` | 30 | overhang angle (°) the deformation aims for. Lower = warps harder, more tilt |
| `ROTATION_MULTIPLIER` | 2 | scales the target tilt. Higher = more tilt |
| `NEIGHBOUR_LOSS_WEIGHT` | 30 | smoothness of the tilt field. Higher = smoother, gentler changes |
| `SET_INITIAL_ROTATION_TO_ZERO` | false | pull areas without overhangs toward no tilt (less noisy) |
| `STEEP_OVERHANG_COMPENSATION` | true | extra rotation for cells that would print in air |
| `INITIAL_ROTATION_FIELD_SMOOTHING` | 30 | 0 = off; any other value = one smoothing pass ([why](../S4_PIPELINE.md#5-known-limitations)) |
| `MAX_POS_ROTATION`, `MAX_NEG_ROTATION` | ±360° | clamp on the target rotation |
| `PART_OFFSET` | [0, 0, 0] | shift the part on the plate, in mm (subtracted) |
| `ROTATION_ITERATIONS` | 100 | solver budget for the tilt field. Changing it changes the result, not just the speed |
| `iterations` | (none) | multi-iteration schedule, see above |

### Turning the field into a shape
| setting | default | what it does |
|---|---|---|
| `DEFORMATION_METHOD` | `island_free` | `island_free` (fit, then lift; [how it works](DEVELOPING.md#the-deformation-island_free-default)) or `notebook` (the notebook's own solve) |
| `DEFORMATION_ITERATIONS` | 1000 | solver budget, `notebook` method only. Changing it changes the result |
| `FIT_METHOD` | `penalty` | `island_free`: `penalty` (robust on any mesh) or `barrier` (strictly fold-free, but can stall on fine meshes) |
| `FOLD_PENALTY`, `FOLD_PENALTY_DET` | 100, 0.5 | `penalty`: how hard tets squashed below `FOLD_PENALTY_DET` × their volume (or folded) are pushed back. Lower lets thin features squash |
| `PENALTY_ITERATIONS` | 300 | `penalty`: solver iterations for the fit |
| `PRECOND_FLOOR` | 0 | solver stiffness floor for the fit (× the median) |
| `ISLAND_LIFT_SLOPE` | 0.5 | every point must be reachable from the bed rising at least this much per mm. Higher = stricter (1.0 ≈ 45° overhangs) but more distortion |
| `ISLAND_LIFT_ROUNDS` | 5 | rounds of lifting |
| `LIFT_ITERATIONS` | 100 | solver iterations per lifting round |
| `LIFT_WEIGHT` | 50 | strength of the lift targets against the fold penalty. Raise it together with `FOLD_PENALTY` |
| `LIFT_ANCHOR` | 1 | while lifting, the grounded vertices stay where the fit put them (0 = off) |
| `LIFT_ONE_SIDED` | true | only the pits are pushed up; earlier lift targets never pull a vertex back down |
| `LIFT_HOLD` | 0.2 | every vertex that isn't in a pit is held where the fit put it, with this fraction of `LIFT_WEIGHT` (0 = off) |
| `LIFT_HOLD_FALLOFF` | 0 | mm over which the hold ramps up with distance from a pit (fewer creases, but costs support on some models) |
| `LIFT_PRECOND_FLOOR` | 0.3 | solver stiffness floor for the lift solves only (`-1` = use `PRECOND_FLOOR`) |
| `BED_TOL` | 0.3 | mm; vertices this close to the lowest point count as the bed face, which is always supported |
| `BED_PIN_WEIGHT` | 0 | > 0 also holds the bed face flat on the bed. Off by default: it costs tilt and support |
| `FLIP_FREE_STAGES`, `FLIP_FREE_STAGE_ITERATIONS`, `BARRIER_WEIGHT` | 10, 150, 0.02 | `barrier` only: tilt ramp stages and barrier strength |
| `SLIVER_QUALITY`, `MICRO_TET_VOLUME` | 0.03, 1e-3 | `barrier` only: badly shaped or tiny tets get no barrier (they would stall it) |

The four `LIFT_ANCHOR` / `LIFT_ONE_SIDED` / `LIFT_HOLD` / `LIFT_PRECOND_FLOOR` switches together are the "anchored,
pit-only" lift. The old lift is `LIFT_ANCHOR=0 LIFT_ONE_SIDED=false LIFT_HOLD=0 LIFT_PRECOND_FLOOR=-1`. Why the new
one is the default, with numbers, is in the
[CHANGELOG](../CHANGELOG.md#2026-10-01-smoother-lifted-surfaces-six-model-quality-suite).

## `map`: the machine and the final 4-axis moves

Notebook cells 17-18, plus the machine limits.

### Geometry and path
| setting | default | what it does |
|---|---|---|
| `NOZZLE_OFFSET` | 42 | mm from the B pivot to the nozzle tip. 42 is the value in the notebook's code; its comment says 41.5. Measure yours on the printer and set it ([the sim has its own copy](../sim/README.md#on-the-real-printer)) |
| `MIN_ROTATION`, `MAX_ROTATION` | −130, 30 | B-axis limits in degrees |
| `ROTATION_AVERAGING_ALPHA` | 0.2 | smoothing of B along the path. Lower = smoother, slower to react |
| `ROTATION_MAX_DELTA` | 1° | B steps bigger than this are split into smaller moves |
| `SEG_SIZE` | 0.6 | mm; long moves are split into segments this long |
| `RETRACTION_LENGTH` | null | mm; `null` = use Cura's `retraction_amount` |

### Extrusion and travel fixes
Each is on by default and off with `--notebook-exact` (see [Differences](#differences-from-the-notebook)).

| setting | default | what it does |
|---|---|---|
| `MAX_EXTRUSION_MULTIPLIER` | 10 | hard cap on the extrusion compensation |
| `EXTRUSION_MULTIPLIER_RANGE` | [0.5, 2.0] | clamp on the extrusion compensation (`null` = off) |
| `SMOOTH_EXTRUSION_MULTIPLIER` | true | smooth the extrusion compensation along the path |
| `SPLIT_RETRACTIONS` | true | retract/unretract in place |
| `SAFE_TRAVEL_TRANSITIONS` | true | lower the nozzle after a travel that left the part, before printing |

### Machine limits
| setting | default | what it does |
|---|---|---|
| `LIMIT_AXIS_SPEEDS` | true | slow down any move that would drive an axis past the limits below |
| `MAX_SPEED_C`, `MAX_SPEED_B` | 360, 180 | deg/s, bed rotation and nozzle tilt. Placeholders: set your machine's real limits |
| `MAX_SPEED_X`, `MAX_SPEED_Z` | 150, 50 | mm/s. Placeholders, as above |
| `MAX_SPEED_E`, `MAX_ACCEL_X`, `MAX_ACCEL_Z`, `MAX_ACCEL_B`, `MAX_ACCEL_C`, `MAX_ACCEL_E`, `CORNER_SPEED` | 60, 2000, 200, 5000, 3000, 3000, 8 | for the print-time estimate only (mm/s, mm/s², deg/s²; `CORNER_SPEED` is the sim's corner speed floor). Placeholders equal to R-Theta Sim's defaults |
| `HOME_X`, `HOME_Z`, `HOME_B`, `HOME_SPEED` | 140, 200, 0, 40 | where `G28` goes and how fast, for the estimate. Same as the sim |

## `cura`: slicing

| setting | default | what it does |
|---|---|---|
| `config` | `cura_config.3mf` | Cura project with the base profile (save one from Cura to use your own) |
| `set` | {} | Cura setting overrides by **internal** name, e.g. `layer_height`, `infill_sparse_density`, `infill_pattern`, `wall_line_count`, `speed_print`, `material_print_temperature` |
| `strip_start_prime` | true | drop Cura's start-code prime |

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

## Differences from the notebook

This is also the answer to "which setting turns on which fix". All of these are on by default. `--notebook-exact`
turns every one of them off and reproduces the notebook's output byte for byte (everything after the
[settings block](../S4_PIPELINE.md#the-settings-block-at-the-top-of-the-g-code)). The rotation (tilt) field is the
notebook's in every case.

| what | setting | difference |
|---|---|---|
| deformation | `DEFORMATION_METHOD` | fold-free, island-free shape instead of the notebook's unconverged solve |
| retractions | `SPLIT_RETRACTIONS` | retract/unretract in place, at Cura's retraction speed (the notebook extruded during a 1 mm plunge, leaving "sticks") |
| extrusion compensation | `SMOOTH_EXTRUSION_MULTIPLIER`, `EXTRUSION_MULTIPLIER_RANGE` | blended smoothly along the path and clamped to 0.5×-2× (was constant per tet, so flow jumped) |
| start code | `strip_start_prime` | Cura's prime is dropped (it became a floating blob inside the part) |
| travel re-entry | `SAFE_TRAVEL_TRANSITIONS` | after a travel that left the part, the nozzle lowers before printing (the notebook drew "poles"); rotation-split steps keep their own move's command; points mapped just below the bed are clamped to it instead of dropped |
| nozzle offset | `NOZZLE_OFFSET` | taken from the profile (was hard-coded) |
| axis speeds | `LIMIT_AXIS_SPEEDS`, `MAX_SPEED_*` | each move takes at least as long as every axis needs at its limit (the notebook used the planar move's time); the G94 `F20000` moves become G93 too |

Why each of these was needed, with measurements, is in the [CHANGELOG](../CHANGELOG.md).
