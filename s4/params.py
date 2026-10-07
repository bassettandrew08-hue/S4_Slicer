"""
Deformation parameters. The notebook keys (PART_OFFSET .. STEEP_OVERHANG_COMPENSATION) default to main.ipynb's
values (cells 2, 4, 7), the ones that produced the verified benchy output. The rest choose and tune the deformation
that replaces the notebook's least-squares solve (DEFORMATION_METHOD, s4/island_free.py). This is the single source
of the deform defaults: island_free.DEFAULTS is built from it.

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
    # Mesh preparation (not in the notebook): coarsen needlessly fine surface tessellation before TetGen
    # (meshio_s4.simplify_surface). 0 = off, the notebook's mesh. ~1e-5 moves the surface ~0.01 mm
    SURFACE_SIMPLIFY_ERROR=0.0,
    # How the rotation field becomes a deformed mesh:
    #   "island_free" (default): fold-free fit + lifting of height minima that would print as floating islands
    #   "notebook": the notebook's least-squares solve (reproduces main.ipynb; folds and leaves islands)
    DEFORMATION_METHOD="island_free",
    # ---- island_free (s4/island_free.py; its module docstring explains the fit and the lift)
    ISLAND_LIFT_SLOPE=0.5,         # mm of rise per mm from the bed; 0.5 ~ overhangs up to ~63 deg from vertical
                                   # in deformed space
    ISLAND_LIFT_ROUNDS=5,
    LIFT_WEIGHT=50.0,
    # Lifting holds the part's grounded vertices (anchor) and everything that is not in a pit (hold) where the fit
    # put them, and only pushes the pits up (one-sided). The old, unanchored lift is LIFT_ANCHOR=0,
    # LIFT_ONE_SIDED=false, LIFT_HOLD=0, LIFT_PRECOND_FLOOR=-1.
    LIFT_ANCHOR=1.0,               # > 0: while lifting, hold the grounded vertices at their height
                                   # (weight x LIFT_WEIGHT)
    LIFT_ONE_SIDED=True,           # lift targets only push up, never hold a vertex down
    LIFT_HOLD=0.2,                 # > 0: hold every vertex that was never in a pit at its fit-only height with a
                                   # two-sided spring (weight x LIFT_WEIGHT), so a pit's rim cannot be dragged up
    LIFT_HOLD_FALLOFF=0.0,         # mm (path length through the mesh): the hold ramps from 0 at the pit to full here
    LIFT_PRECOND_FLOOR=0.3,        # PRECOND_FLOOR for the lifting solves only (< 0: same as the fit's). Lift
                                   # springs pull the whole mesh; without a floor the preconditioned step is huge at
                                   # vertices held only by micro-tets, and a dense micro-tet cluster gets torn apart
    BARRIER_WEIGHT=0.02,           # barrier: beta, weight of the symmetric-Dirichlet fold barrier
    FLIP_FREE_STAGES=10,           # barrier: the target rotations are ramped up in this many stages
    FLIP_FREE_STAGE_ITERATIONS=150,
    FIT_METHOD="penalty",          # "penalty": exact linear fit, then a soft fold penalty (robust on any mesh);
                                   # "barrier": fold-free barrier + staged homotopy (strict, but can stall on fine
                                   # meshes)
    FOLD_PENALTY=100.0,            # penalty: strength of the soft anti-fold term
    FOLD_PENALTY_DET=0.5,          # penalty: tets squashed below this volume ratio (or folded) are penalised
    PENALTY_ITERATIONS=300,
    LIFT_ITERATIONS=100,           # solver iterations per lifting round (warm-started)
    BED_PIN_WEIGHT=0.0,            # > 0: hold the part's bed face flat on the bed. Off: pinning fights the tilt near
                                   # the base (the benchy got 1475 mm of unsupported extrusion at 500)
    BED_TOL=0.3,                   # mm: vertices this close to the lowest point form the bed face (always grounded)
    SLIVER_QUALITY=0.03,           # barrier: tets with mean-ratio shape quality below this get no fold barrier
                                   # (they would stall the solve)
    MICRO_TET_VOLUME=1e-3,         # barrier: ... nor do tets smaller than this fraction of the median tet volume
    PRECOND_FLOOR=0.0,             # minimum per-vertex stiffness in the solver's preconditioner (x median)
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
