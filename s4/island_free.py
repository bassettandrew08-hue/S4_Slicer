"""
Island-free, fold-free deformation (DEFORMATION_METHOD = "island_free").

The notebook turns the optimised rotation field into a deformed mesh with an unconverged least-squares solve.
That solve folds the mesh (202 inverted tets on the benchy) and, more importantly, leaves local height minima
in the deformed shape. A planar slicer starts printing every such minimum in mid-air: those are the floating
islands. This module replaces only that step; the rotation field is still the notebook's (unless PIT_FREE_ROUNDS).

1. Fold-free fit (FIT_METHOD "barrier"; the default "penalty" replaces the barrier below by a soft penalty on
   tets squashed below FOLD_PENALTY_DET, see PenaltyProblem). Minimise over vertex positions V
       sum_c vol_c ( |J_c - R_c|^2  +  beta (|J_c|^2 + |J_c^-1|^2 - 6) )
   J_c = deformation gradient of tet c, R_c = its target rotation. The second term (symmetric Dirichlet) is a
   barrier that is infinite at zero volume, and every step is capped before the first tet would flip, so the mesh
   can never fold. Large rotations are reached by ramping the targets up in stages (homotopy). L-BFGS,
   preconditioned with the mesh stiffness matrix.
2. Island lifting. Priority-flood from the bed gives, for every vertex, the lowest height at which it can be
   reached from the bed by a path rising at least ISLAND_LIFT_SLOPE per mm. Vertices below that height sit in a
   "pit" that would print in mid-air; they get soft height targets and the fit is re-solved, a few rounds. The
   targets only push up (LIFT_ONE_SIDED), the grounded vertices are anchored (LIFT_ANCHOR), and every vertex that
   isn't in a pit is held at its fit-only height (LIFT_HOLD), so only the pits move: an unanchored lift mostly moved
   the whole part up, and its stale targets pulled vertices back down (creases). The lift solves use their own
   preconditioner floor (LIFT_PRECOND_FLOOR) so vertices held only by micro-tets don't take huge steps.
   LIFT_MODE "band" (opt-in) fills the pits differently: the flood level is a cone around the pit's outlet, so on a
   wall the layer lines close into rings around the lifted patch. "band" instead gives the pit (and a margin of
   LIFT_BAND_MARGIN mm) two-sided targets whose surface heights are the harmonic interpolation of the heights on the
   patch's rim (band_fill): no new extrema on the surface, so the layer lines run through the patch as bands. It
   is meant for a low ISLAND_LIFT_SLOPE (~0.15): a steep slope forces the cone back in.

Opt-in (PIT_FREE_ROUNDS > 0), between 1 and 2: pit_free_field changes the rotation field itself where the fit
leaves pits, until there are none (the one place the notebook's field is changed; off by default).

The solvers (numba kernels, FitProblem, PenaltyProblem) are in s4/island_free_solver.py. Deterministic: all
per-tet work is written to arrays and reduced serially.
"""
import heapq
import time

import numpy as np
from scipy.spatial.transform import Rotation

from .geometry import rotation_matrices, tangential_vectors
from .island_free_solver import FitProblem, PenaltyProblem, mean_ratio  # noqa: F401 (re-exported)
from .params import DEFAULT_PARAMS

# the deform settings this module reads; defaults, order and their explanations are in s4/params.py
_KEYS = (
    "ISLAND_LIFT_SLOPE", "ISLAND_LIFT_ROUNDS", "LIFT_MODE", "LIFT_BAND_MARGIN", "LIFT_WEIGHT", "LIFT_ANCHOR", "LIFT_ONE_SIDED", "LIFT_HOLD",
    "LIFT_HOLD_FALLOFF", "LIFT_PRECOND_FLOOR", "BARRIER_WEIGHT", "FLIP_FREE_STAGES", "FLIP_FREE_STAGE_ITERATIONS",
    "FIT_METHOD", "FOLD_PENALTY", "FOLD_PENALTY_DET", "PENALTY_ITERATIONS", "LIFT_ITERATIONS", "BED_PIN_WEIGHT",
    "BED_TOL", "SLIVER_QUALITY", "MICRO_TET_VOLUME", "PRECOND_FLOOR",
    "PIT_FREE_ROUNDS", "PIT_FREE_SLOPE", "PIT_FREE_MAX_STEP", "PIT_FREE_SMOOTH", "PIT_FREE_FRACTION",
    "PIT_FREE_REGION", "PIT_FREE_HOLD", "PIT_FREE_SEED_DEPTH",
)
DEFAULTS = {k: DEFAULT_PARAMS[k] for k in _KEYS}


# ---------------------------------------------------------------------------------------- island tools

def vertex_graph(cells, n):
    e = np.vstack([cells[:, [i, j]] for i in range(4) for j in range(i + 1, 4)])
    e = np.unique(np.sort(e, axis=1), axis=0)
    nbrs = [[] for _ in range(n)]
    for a, b in e.tolist():
        nbrs[a].append(b)
        nbrs[b].append(a)
    return nbrs


def priority_flood(V, nbrs, slope, bed_tol, bed=None, spill=False):
    """Lowest heights >= z so that every vertex is reachable from the bed by a path rising >= slope per mm. The bed
    is the given vertices (the pinned bed face) or else everything within bed_tol of the lowest point.
    spill=True also returns, per vertex, the vertex whose own height set its level (itself unless it is in a pit;
    for a pit vertex the rim vertex the pit spills over)."""
    z = V[:, 2]
    n = len(z)
    h = np.full(n, np.inf)
    src = np.arange(n)
    done = np.zeros(n, bool)
    pq = []
    for v in (np.nonzero(z <= z.min() + bed_tol)[0] if bed is None else bed).tolist():
        h[v] = z[v]
        heapq.heappush(pq, (float(z[v]), v))
    xy = V[:, :2].tolist()
    zl = z.tolist()
    while pq:
        hv, v = heapq.heappop(pq)
        if done[v]:
            continue
        done[v] = True
        xv, yv = xy[v]
        for u in nbrs[v]:
            if done[u]:
                continue
            xu, yu = xy[u]
            rise = hv + slope * ((xu - xv) ** 2 + (yu - yv) ** 2) ** 0.5
            cand = max(zl[u], rise)
            if cand < h[u]:
                h[u] = cand
                src[u] = u if zl[u] >= rise else src[v]
                heapq.heappush(pq, (cand, u))
    return (h, src) if spill else h


def island_seeds(V, nbrs, bed_tol=1.0, min_persistence=0.2, bed=None):
    """Local height minima not connected to the bed, with how much height they float for (sublevel-set persistence).
    Returns [(birth vertex, birth z, merge z, vertices)]."""
    z = V[:, 2]
    n = len(z)
    parent = list(range(n))
    birth = z.copy()
    size = np.ones(n, int)
    if bed is None:
        grounded = z <= z.min() + bed_tol
    else:
        grounded = np.zeros(n, bool)
        grounded[bed] = True

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    active = np.zeros(n, bool)
    seeds = []
    for v in np.argsort(z, kind="stable").tolist():
        active[v] = True
        for u in nbrs[v]:
            if not active[u]:
                continue
            ru, rv = find(u), find(v)
            if ru == rv:
                continue
            if grounded[ru] and not grounded[rv]:
                young, old = rv, ru
            elif grounded[rv] and not grounded[ru]:
                young, old = ru, rv
            else:
                young, old = (rv, ru) if birth[rv] > birth[ru] else (ru, rv)
            if not grounded[young] and z[v] - birth[young] >= min_persistence:
                seeds.append((int(young), float(birth[young]), float(z[v]), int(size[young])))
            parent[young] = old
            size[old] += size[young]
            grounded[old] = grounded[old] or grounded[young]
            birth[old] = min(birth[old], birth[young])
    return seeds


# ----------------------------------------------------------------------------------------------- driver

def _grounded(V, bed, tol=0.5):
    """Vertices standing on the bed: the part's bed face, plus anything within tol of the lowest point."""
    return np.union1d(bed, np.nonzero(V[:, 2] <= V[:, 2].min() + tol)[0]).astype(np.int64)


def _hold_targets(V, nbrs, pit, held, weight, falloff):
    """Two-sided z springs (per-vertex weight) at V's height for the held vertices outside the pit, ramping
    (smoothstep) from 0 at the pit to full weight at a path length of falloff mm through the mesh."""
    keep = held.copy()
    keep[pit] = False
    idx = np.nonzero(keep)[0].astype(np.int64)
    w = weight[idx].copy()
    if falloff > 0 and len(pit):
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import dijkstra
        a = np.repeat(np.arange(len(nbrs)), [len(x) for x in nbrs])
        b = np.concatenate([np.asarray(x, np.int64) for x in nbrs])
        G = csr_matrix((np.linalg.norm(V[a] - V[b], axis=1) + 1e-9, (a, b)), shape=(len(nbrs), len(nbrs)))
        d = dijkstra(G, directed=False, indices=pit, min_only=True, limit=falloff)
        t = np.clip(d[idx] / falloff, 0.0, 1.0)
        w *= t * t * (3.0 - 2.0 * t)
    return idx, V[idx, 2].copy(), w, np.zeros(len(idx), bool)


def pit_free_field(prob, cell_centers, rotation_field, V, bed, nbrs, q, solve):
    """Pit-aware rotation field (PIT_FREE_ROUNDS > 0): change the rotation field where the fit leaves pits.

    The fit's deformed height is (to first order, the linear fit) h = K^-1 b(theta), b = sum_c vol_c G_c g_c with
    g_c = cos(theta_c) z + sin(theta_c) (z x t_c), the third row of the tet's target rotation (G_c: its basis-function
    gradients, K the stiffness matrix, t_c its tangential axis). Each round:
      1. the priority flood (slope PIT_FREE_SLOPE) gives every pit vertex its deficit d_v (how far it sits below the
         level at which it is reachable from the bed) and the rim vertex its pit spills over;
      2. the pit depth P = sum_v d_v (h_rim(v) - h_v) is differentiated with respect to every tet's angle with one
         adjoint solve, lam = K^-1 a (a = -d_v at the pit vertices, +d_v at their rims):
         dP/dtheta_c = vol_c sum_k lam_ck (G_ck . dg_c/dtheta_c);
      3. that gradient (per unit volume) is spread smoothly, (M + s^2 K) u = M g on the vertices (M lumped volume,
         s = PIT_FREE_SMOOTH mm), back to the tets as the mean of their vertices;
      4. theta -= a u, a chosen so the linear model removes PIT_FREE_FRACTION of P, but no tet's angle changes by more
         than PIT_FREE_MAX_STEP deg; then the fit is solved again for the new field (solve(R) -> V).
    With PIT_FREE_REGION > 0 the change is confined to a smooth mask around the first fit's pits and the vertices
    beyond it are held at their heights (solve(R, hold)), so the rest of the part keeps its shape.
    Stops when no vertex is more than 0.02 mm deep in a pit. Returns (rotation field, R, V, info)."""
    from scipy.sparse import diags, identity
    from scipy.sparse.linalg import splu
    n = prob.n
    cells = prob.cells
    Gi = prob.Dm_inv
    G = np.concatenate([-Gi.sum(axis=1, keepdims=True), Gi], axis=1)
    K = prob.K_raw
    slope = float(q["PIT_FREE_SLOPE"])
    region_mm = float(q["PIT_FREE_REGION"])
    mask = hold = None
    wz = np.zeros(n)
    if region_mm > 0:
        # local change: only the tets within PIT_FREE_REGION mm (path length through the mesh) of the fit's pits may
        # change, full weight up to half of that and smoothly down to 0 at the edge (a tet takes the smallest weight
        # of its vertices); every vertex beyond the edge is held at its height (z spring, PIT_FREE_HOLD), so the
        # rest of the part keeps the fit's shape. The adjoint uses the same held stiffness.
        from scipy.sparse.csgraph import dijkstra
        h0 = priority_flood(V, nbrs, slope, bed_tol=0.5, bed=_grounded(V, bed))
        _check_flood(h0, "PIT_FREE")
        pit0 = np.nonzero(h0 - V[:, 2] > max(0.02, float(q["PIT_FREE_SEED_DEPTH"])))[0]
        if not len(pit0):
            # no pit deeper than PIT_FREE_SEED_DEPTH: nothing to confine the change to, so the field stays as it is
            # (without a mask the rounds would change the whole field)
            d0 = h0 - V[:, 2]
            p0 = np.nonzero(d0 > 0.02)[0]
            depth0 = round(float(d0[p0].sum()), 2)
            info = {"pit_vertices_before": len(p0), "pit_depth_before": depth0, "pit_vertices": len(p0),
                    "pit_depth": depth0, "pit_free_rounds": 0}
            rf = np.asarray(rotation_field, np.float64).copy()
            return rf, rotation_matrices(cell_centers, rf), V, info
        if len(pit0):
            dist = dijkstra(_edge_graph(prob.P0, nbrs), directed=False, indices=pit0, min_only=True, limit=region_mm)
            t = np.clip((region_mm - dist) / (0.5 * region_mm), 0.0, 1.0)
            mask = (t * t * (3.0 - 2.0 * t))[cells].min(axis=1)
            out = np.nonzero(~np.isfinite(dist))[0].astype(np.int64)
            if float(q["PIT_FREE_HOLD"]) > 0:
                hold = (out, V[out, 2].copy(), np.full(len(out), float(q["PIT_FREE_HOLD"])))
                wz[out] = float(q["PIT_FREE_HOLD"])
    lu = splu((K + diags(wz + 1e-9 * np.median(K.diagonal()))).tocsc())
    zt = np.cross(np.array([0.0, 0.0, 1.0]), tangential_vectors(cell_centers))
    s = float(q["PIT_FREE_SMOOTH"])
    smooth = None
    if s > 0:
        Mv = np.zeros(n)
        np.add.at(Mv, cells.ravel(), np.repeat(prob.vol / 4, 4))
        smooth = splu((diags(Mv) + s * s * K).tocsc())
    rf = np.asarray(rotation_field, np.float64).copy()
    rounds = int(q["PIT_FREE_ROUNDS"])
    R = rotation_matrices(cell_centers, rf)
    first = None
    k = 0
    for k in range(rounds + 1):
        h, src = priority_flood(V, nbrs, slope, bed_tol=0.5, bed=_grounded(V, bed), spill=True)
        _check_flood(h, "PIT_FREE")
        d = h - V[:, 2]
        pit = np.nonzero(d > 0.02)[0]
        depth = float(d[pit].sum())
        if first is None:
            first = (len(pit), depth)
        if not len(pit) or k == rounds:
            break
        a = np.zeros(n)
        np.add.at(a, pit, -d[pit])
        np.add.at(a, src[pit], d[pit])
        lam = lu.solve(a)
        dg = -np.sin(rf)[:, None] * np.array([0.0, 0.0, 1.0]) + np.cos(rf)[:, None] * zt
        grad = prob.vol * np.einsum("ck,ck->c", np.einsum("ckj,cj->ck", G, dg), lam[cells])  # dP/dtheta
        u = grad / prob.vol
        if smooth is not None:
            gv = np.zeros(n)
            np.add.at(gv, cells.ravel(), np.repeat(grad / 4, 4))
            u = smooth.solve(gv)[cells].mean(axis=1)
        if mask is not None:
            u = u * mask
        umax = float(np.abs(u).max())
        gu = float(grad @ u)
        if umax <= 0 or gu <= 0:
            break
        step = min(float(q["PIT_FREE_FRACTION"]) * depth / gu, np.radians(float(q["PIT_FREE_MAX_STEP"])) / umax)
        rf -= step * u
        R = rotation_matrices(cell_centers, rf)
        V = solve(R, hold)
    info = {"pit_vertices_before": first[0], "pit_depth_before": round(first[1], 2), "pit_vertices": len(pit),
            "pit_depth": round(depth, 2), "pit_free_rounds": k}
    return rf, R, V, info


def _check_flood(h, what):
    """The priority flood reaches every vertex connected to the bed through the mesh; an infinite level means a body
    that is not (a separate shell floating above the bed), which no height target can fix."""
    bad = ~np.isfinite(h)
    if bad.any():
        raise ValueError(f"{what}: {int(bad.sum())} mesh vertices are not connected to the bed through the mesh "
                         f"(a separate body that does not touch the bed?); the island lift cannot place them. "
                         f"Join the bodies or move them onto the bed.")


def _edge_graph(P, nbrs):
    """Sparse symmetric graph of the mesh edges weighted by their length in P."""
    from scipy.sparse import csr_matrix
    a = np.repeat(np.arange(len(nbrs)), [len(x) for x in nbrs])
    b = np.concatenate([np.asarray(x, np.int64) for x in nbrs])
    return csr_matrix((np.linalg.norm(P[a] - P[b], axis=1) + 1e-9, (a, b)), shape=(len(nbrs), len(nbrs)))


def surface_edges(cells):
    """Edges of the tet mesh's boundary triangles."""
    c = np.asarray(cells)
    fs = np.sort(np.vstack([c[:, [0, 1, 2]], c[:, [0, 1, 3]], c[:, [0, 2, 3]], c[:, [1, 2, 3]]]), axis=1)
    u, cnt = np.unique(fs, axis=0, return_counts=True)
    sf = u[cnt == 1]
    return np.unique(np.sort(np.vstack([sf[:, [0, 1]], sf[:, [1, 2]], sf[:, [0, 2]]]), axis=1), axis=0)


def harmonic_fill(z, region, P, edges, fallback=None):
    """z with the region's values replaced by the harmonic interpolation (1/length^2 edge weights over `edges`, which
    obey the discrete maximum principle: no minima or maxima inside) of z on the region's boundary. A connected part
    of the region with no boundary (a whole surface shell inside it) has no rim to interpolate from (its Laplacian is
    singular): it keeps `fallback` (default z)."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse.linalg import spsolve
    n = len(z)
    a, b = edges[:, 0], edges[:, 1]
    keep = region[a] | region[b]
    a, b = a[keep], b[keep]
    w = 1.0 / np.maximum(np.linalg.norm(P[a] - P[b], axis=1), 1e-6) ** 2
    W = csr_matrix((np.r_[w, w], (np.r_[a, b], np.r_[b, a])), shape=(n, n))
    D = np.nonzero(region & (np.asarray(W.sum(1)).ravel() > 0))[0]
    f = z.copy()
    if not len(D):
        return f
    Wd = W[D]
    L = csr_matrix((np.asarray(Wd.sum(1)).ravel(), (np.arange(len(D)), np.arange(len(D)))),
                   shape=(len(D), len(D))) - Wd[:, D]
    out = np.ones(n, bool)
    out[D] = False
    rhs = Wd[:, out] @ z[out]
    _, comp = connected_components(Wd[:, D], directed=False)
    rim = np.zeros(comp.max() + 1, bool)
    rim[comp[np.asarray(Wd[:, out].sum(1)).ravel() > 0]] = True
    ok = rim[comp]
    if not ok.all():
        f[D[~ok]] = (z if fallback is None else fallback)[D[~ok]]
        if not ok.any():
            return f
        L = L.tocsr()[ok][:, ok]
        D, rhs = D[ok], rhs[ok]
    f[D] = spsolve(L.tocsc(), rhs)
    return f


def band_fill(P0, V, nbrs, sedges, region, h, slope, bed):
    """LIFT_MODE "band": heights for the vertices of `region` (bool: the pits and a margin around them) whose layer
    lines continue the ones around them instead of closing into rings. On the outer surface: the harmonic
    interpolation, over the surface triangles, of the heights on the region's rim (no new extrema on the surface, so
    no closed rings); inside the part: the flood level h. The result is flooded again at `slope` from the bed, so
    anything still unreachable rises to the flood level. Returns the full height field (V's z outside the region)."""
    z = V[:, 2]
    on_s = np.zeros(len(z), bool)
    on_s[sedges.ravel()] = True
    f = np.where(region, np.maximum(z, h), z)
    f = np.where(region & on_s, harmonic_fill(z, region & on_s, P0, sedges, fallback=f), f)
    Vf = V.copy()
    Vf[:, 2] = f
    return np.maximum(f, priority_flood(Vf, nbrs, slope, bed_tol=0.5, bed=bed))


def deform(points, cells, cell_centers, rotation_field, p, lift=True, log=None):
    """Fold-free fit of the rotation field, then island lifting. Returns (new vertex positions, info dict)."""
    q = {**DEFAULTS, **{k: p[k] for k in DEFAULTS if k in p}}
    t0 = time.perf_counter()
    P0 = np.asarray(points, dtype=np.float64)
    V = P0.copy()
    iters = 0
    # the part's bed face: supported by the bed wherever the deformation puts it
    bed = np.nonzero(P0[:, 2] <= P0[:, 2].min() + float(q["BED_TOL"]))[0].astype(np.int64)
    pin = None
    if float(q["BED_PIN_WEIGHT"]) > 0:
        pin = (bed, P0[bed, 2].copy(), np.full(len(bed), float(q["BED_PIN_WEIGHT"])))
    if q["FIT_METHOD"] == "penalty":
        prob = PenaltyProblem(P0, cells, gamma=float(q["FOLD_PENALTY"]), eps=float(q["FOLD_PENALTY_DET"]),
                              precond_floor=float(q["PRECOND_FLOOR"]))
        R = rotation_matrices(cell_centers, rotation_field)
        V = prob.linear_start(R, pin)
        V, iters = prob.solve(R, V, int(q["PENALTY_ITERATIONS"]), lift=pin)
    elif q["FIT_METHOD"] == "barrier":
        prob = FitProblem(P0, cells, q["BARRIER_WEIGHT"], sliver_quality=float(q["SLIVER_QUALITY"]),
                          micro_volume=float(q["MICRO_TET_VOLUME"]), precond_floor=float(q["PRECOND_FLOOR"]))
        stages = int(q["FLIP_FREE_STAGES"])
        R = None
        for k in range(1, stages + 1):
            R = rotation_matrices(cell_centers, (k / stages) * rotation_field)
            V, it = prob.solve(R, V, int(q["FLIP_FREE_STAGE_ITERATIONS"]), lift=pin)
            iters += it
    else:
        raise ValueError(f"unknown FIT_METHOD {q['FIT_METHOD']!r} (use 'penalty' or 'barrier')")
    t_fit = time.perf_counter() - t0
    nbrs = vertex_graph(np.asarray(cells), len(P0))
    pit_info = None
    rotation_field_in = rotation_field  # Byrd's field (PIT_FREE_ROUNDS replaces rotation_field below)
    if int(q["PIT_FREE_ROUNDS"]) > 0:
        if q["FIT_METHOD"] != "penalty":
            raise ValueError("PIT_FREE_ROUNDS needs FIT_METHOD 'penalty'")
        rf_in = np.asarray(rotation_field, np.float64)

        def _solve(Rk, hold=None):
            nonlocal iters
            tg = pin if hold is None else (hold if pin is None else tuple(np.concatenate([a, b])
                                                                           for a, b in zip(pin, hold)))
            Vk, it = prob.solve(Rk, prob.linear_start(Rk, tg), int(q["PENALTY_ITERATIONS"]), lift=tg)
            iters += it
            return Vk
        rotation_field, R, V, pit_info = pit_free_field(prob, cell_centers, rf_in, V, bed, nbrs, q, _solve)
        dd = np.abs(np.degrees(rotation_field - rf_in))
        pit_info["field_change_deg_mean"] = round(float((dd * prob.vol).sum() / prob.vol.sum()), 2)
        pit_info["field_change_deg_max"] = round(float(dd.max()), 1)
        pit_info["pit_free_seconds"] = round(time.perf_counter() - t0 - t_fit, 2)
    lifted = 0
    if lift:
        tgt = np.full(len(P0), -np.inf)
        w = np.full(len(P0), float(q["LIFT_WEIGHT"]))
        anchor = None
        if float(q["LIFT_ANCHOR"]) > 0:
            g0 = _grounded(V, bed)
            anchor = (g0, V[g0, 2].copy(), float(q["LIFT_ANCHOR"]) * w[g0], np.zeros(len(g0), bool))
        if float(q["LIFT_PRECOND_FLOOR"]) >= 0 and float(q["LIFT_PRECOND_FLOOR"]) != float(q["PRECOND_FLOOR"]):
            prob.set_precond_floor(float(q["LIFT_PRECOND_FLOOR"]))
        # no hold together with the bed pin: both fix z near the base and the lift squeezes the tets between them
        hold_w = np.full(len(P0), 0.0 if pin is not None else float(q["LIFT_HOLD"]) * float(q["LIFT_WEIGHT"]))
        V_fit = V.copy()
        held = np.ones(len(P0), bool)
        if anchor is not None:
            held[anchor[0]] = False
        if pin is not None:
            held[pin[0]] = False
        band = q["LIFT_MODE"] == "band"
        if q["LIFT_MODE"] not in ("flood", "band"):
            raise ValueError(f"unknown LIFT_MODE {q['LIFT_MODE']!r} (use 'flood' or 'band')")
        if band:
            G0 = _edge_graph(P0, nbrs)
            sedges = surface_edges(cells)
            pit_all = np.zeros(len(P0), bool)
        for _ in range(int(q["ISLAND_LIFT_ROUNDS"])):
            g = _grounded(V, bed)
            h = priority_flood(V, nbrs, float(q["ISLAND_LIFT_SLOPE"]), bed_tol=0.5, bed=g)
            _check_flood(h, "island lift")
            need = (h - V[:, 2]) > 0.02
            if not need.any():
                break
            if band:
                from scipy.sparse.csgraph import dijkstra
                pit_all |= need
                region = pit_all.copy()
                if float(q["LIFT_BAND_MARGIN"]) > 0:
                    d = dijkstra(G0, directed=False, indices=np.nonzero(pit_all)[0], min_only=True,
                                 limit=float(q["LIFT_BAND_MARGIN"]))
                    region |= np.isfinite(d)
                if anchor is not None:
                    region[anchor[0]] = False
                f = band_fill(P0, V, nbrs, sedges, region, h, float(q["ISLAND_LIFT_SLOPE"]), g)
                two = region.copy()
                up = ~region & (f > V[:, 2] + 0.02)
                idx = np.nonzero(two | up)[0].astype(np.int64)
                targets = (idx, f[idx], w[idx], up[idx] & bool(q["LIFT_ONE_SIDED"]))
            else:
                tgt = np.maximum(tgt, np.where(need, h, -np.inf))
                idx = np.nonzero(np.isfinite(tgt))[0].astype(np.int64)
                targets = (idx, tgt[idx], w[idx], np.full(len(idx), bool(q["LIFT_ONE_SIDED"])))
            hold = None
            if float(q["LIFT_HOLD"]) > 0:
                hold = _hold_targets(V_fit, nbrs, idx, held, hold_w, float(q["LIFT_HOLD_FALLOFF"]))
            for extra in (anchor, pin, hold):
                if extra is not None:
                    if len(extra) == 3:
                        extra = extra + (np.zeros(len(extra[0]), bool),)
                    targets = tuple(np.concatenate([a, b]) for a, b in zip(targets, extra))
            V, it = prob.solve(R, V, int(q["LIFT_ITERATIONS"]) if q["FIT_METHOD"] == "penalty" else 300,
                               lift=targets)
            iters += it
            lifted = len(idx)
    seeds = island_seeds(V, nbrs, bed=_grounded(V, bed, tol=1.0))
    c = np.asarray(cells)
    J = np.stack([V[c[:, k]] - V[c[:, 0]] for k in (1, 2, 3)], axis=2) @ prob.Dm_inv
    det = np.linalg.det(J)
    U, _, Wt = np.linalg.svd(J)
    flip = np.linalg.det(U @ Wt) < 0
    U[flip, :, 2] *= -1
    rotvec = Rotation.from_matrix(U @ Wt).as_rotvec()
    about_t = np.degrees(np.einsum("ij,ij->i", rotvec, tangential_vectors(cell_centers)))
    want = np.degrees(rotation_field)
    sel = np.abs(want) > 20.0
    def _wmed(x, w):
        o = np.argsort(x)
        cw = np.cumsum(w[o])
        return float(x[o][np.searchsorted(cw, cw[-1] / 2)])
    # tilt that matters for the B axis: rotation about the tangential axis, volume-weighted, where >20 deg is wanted
    achieved = abs(_wmed(about_t[sel], prob.vol[sel])) if sel.any() else 0.0
    wanted = abs(_wmed(want[sel], prob.vol[sel])) if sel.any() else 0.0  # the field the fit aimed for
    want0 = np.degrees(np.asarray(rotation_field_in, np.float64))
    sel0 = np.abs(want0) > 20.0
    wanted0 = abs(_wmed(want0[sel0], prob.vol[sel0])) if sel0.any() else 0.0  # Byrd's field
    info = {"seconds": round(time.perf_counter() - t0, 2), "fit_seconds": round(t_fit, 2), "iterations": iters,
            "inverted_tets": int((det <= 0).sum()),
            "inverted_volume_pct": round(100.0 * float(prob.vol[det <= 0].sum() / prob.vol.sum()), 4),
            "tilt_deg": round(achieved, 1), "target_tilt_deg": round(wanted0, 1),
            "lifted_vertices": lifted, "island_seeds": len(seeds),
            "island_mass": round(float(sum((d - b) * m for _, b, d, m in seeds)), 1)}
    if pit_info is not None:
        info.update(pit_info)
        info["pit_free_target_tilt_deg"] = round(wanted, 1)  # the changed field's; target_tilt_deg is Byrd's
    if log:
        log(f"[deform] island-free: {info}")
        if wanted > 5 and achieved < 0.7 * wanted:
            log(f"[deform] WARNING: the deformation reached only {achieved:.0f} deg of the {wanted:.0f} deg tilt it "
                f"aimed for (volume-weighted, where >20 deg is wanted). The nozzle will tilt too little; "
                f"try FIT_METHOD \"penalty\".")
    return V, info
