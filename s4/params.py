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
    LIFT_WEIGHT=5.0,
    BARRIER_WEIGHT=0.02,
    FLIP_FREE_STAGES=10,
    FLIP_FREE_STAGE_ITERATIONS=150,
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
