"""
Deformation parameters. Defaults are the main.ipynb values that produced the verified benchy output
(cells 2, 4, 7 as of 2026-09-25 19:00).

These are the "deform" section of a build profile (s4/profile.py). An "iterations" list replays the notebook's
multi-iteration workflow (cell 4 -> cell 7 -> cell 9 "undeformed_tet = deformed_tet.copy()" -> cell 4 ...). Like
notebook variables, each iteration inherits the previous iteration's values and overrides only the keys it lists.
"""
import numpy as np

DEFAULT_PARAMS = dict(
    PART_OFFSET=[0., 0., 0.],
    NEIGHBOUR_LOSS_WEIGHT=30,  # the larger the weight, the more the rotation field will be smoothed
    MAX_OVERHANG=30,  # the maximum overhang angle in degrees
    ROTATION_MULTIPLIER=2,  # the larger the multiplier, the more the rotation field will be rotated
    SET_INITIAL_ROTATION_TO_ZERO=False,
    INITIAL_ROTATION_FIELD_SMOOTHING=30,
    MAX_POS_ROTATION=float(np.deg2rad(360)),  # user fix: was 3600
    MAX_NEG_ROTATION=float(np.deg2rad(-360)),  # user fix: was -3600
    ROTATION_ITERATIONS=100,
    DEFORMATION_ITERATIONS=1000,
    STEEP_OVERHANG_COMPENSATION=True,
    # How the rotation field becomes a deformed mesh:
    #   "island_free" (default): fold-free fit + lifting of height minima that would print as floating islands
    #   "notebook": the notebook's least-squares solve (reproduces main.ipynb; folds and leaves islands)
    DEFORMATION_METHOD="island_free",
    ISLAND_LIFT_SLOPE=0.5,         # island_free: minimum rise per mm from the bed (0.5 ~ 63 deg overhang)
    ISLAND_LIFT_ROUNDS=5,
    LIFT_WEIGHT=50.0,
    LIFT_ANCHOR=1.0,               # island_free: > 0 holds the grounded vertices still while lifting (x LIFT_WEIGHT)
    LIFT_ONE_SIDED=True,           # island_free: lift targets only push up, never hold a vertex down
    LIFT_HOLD=0.2,                 # island_free: > 0 holds never-in-a-pit vertices at their fit height (x LIFT_WEIGHT)
    LIFT_HOLD_FALLOFF=0.0,         # island_free: mm over which that hold ramps up with distance from the pit
    LIFT_PRECOND_FLOOR=0.3,        # island_free: PRECOND_FLOOR for the lifting solves only (< 0 = the fit's)
    BARRIER_WEIGHT=0.02,
    FLIP_FREE_STAGES=10,
    FLIP_FREE_STAGE_ITERATIONS=150,
    FIT_METHOD="penalty",          # island_free: "penalty" (robust, default) or "barrier" (strictly fold-free)
    FOLD_PENALTY=100.0,            # penalty: strength of the soft anti-fold term
    FOLD_PENALTY_DET=0.5,          # penalty: tets squashed below this volume ratio are penalised
    PENALTY_ITERATIONS=300,
    LIFT_ITERATIONS=100,           # island_free: solver iterations per lifting round
    BED_PIN_WEIGHT=0.0,            # island_free: > 0 holds the bed face flat on the bed (off: costs tilt and support)
    BED_TOL=0.3,                   # island_free: mm above the lowest point that counts as the bed face (always supported)
    SLIVER_QUALITY=0.03,           # barrier: badly shaped tets get no fold barrier (they would stall the solve)
    MICRO_TET_VOLUME=1e-3,         # barrier: nor do tets below this fraction of the median volume
    PRECOND_FLOOR=0.0,             # solver preconditioner stiffness floor (x median)
)


def expand_iterations(params=None):
    """-> (part_offset, [per-iteration param dicts])."""
    params = dict(params or {})
    iters = params.pop("iterations", None)
    base = dict(DEFAULT_PARAMS)
    base.update(params)
    if not iters:
        return base["PART_OFFSET"], [base]
    out = []
    cur = base
    for it in iters:
        cur = {**cur, **it}
        out.append(cur)
    return base["PART_OFFSET"], out
