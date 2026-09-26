"""
Deformation parameters. Defaults are the main.ipynb values that produced the verified benchy output
(cells 2, 4, 7 as of 2026-09-25 19:00).

A params JSON may be a flat dict (one deformation iteration), or
    {"PART_OFFSET": [...], "iterations": [{...}, {...}, ...]}
for the notebook's multi-iteration workflow (cell 4 -> cell 7 -> cell 9 "undeformed_tet = deformed_tet.copy()"
-> cell 4 ...). Like notebook variables, each iteration inherits the previous iteration's values and overrides
only the keys it lists.
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
        for k in ("MAX_POS_ROTATION", "MAX_NEG_ROTATION"):
            if isinstance(cur[k], str) and cur[k].startswith("deg:"):
                cur[k] = float(np.deg2rad(float(cur[k][4:])))
        out.append(cur)
    return base["PART_OFFSET"], out
