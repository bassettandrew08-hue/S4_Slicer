"""
Island-free, fold-free deformation (DEFORMATION_METHOD = "island_free").

The notebook turns the optimised rotation field into a deformed mesh with an unconverged least-squares solve.
That solve folds the mesh (202 inverted tets on the benchy) and, more importantly, leaves local height minima
in the deformed shape. A planar slicer starts printing every such minimum in mid-air: those are the floating
islands. This module replaces only that step; the rotation field is still the notebook's.

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

Deterministic: all per-tet work is written to arrays and reduced serially.
"""
import heapq
import time

import numba
import numpy as np
from scipy.sparse import coo_matrix, diags, identity
from scipy.sparse.linalg import splu
from scipy.spatial.transform import Rotation

from .fast_map import rotation_matrices, tangential_vectors

DEFAULTS = dict(
    ISLAND_LIFT_SLOPE=0.5,       # mm of rise per mm; 0.5 ~ overhangs up to ~63 deg from vertical in deformed space
    ISLAND_LIFT_ROUNDS=5,
    LIFT_WEIGHT=50.0,
    # Lifting: hold the part's grounded vertices (anchor) and everything that is not in a pit (hold) where the fit
    # put them, and only push the pits up (one-sided). Without the anchor the lift mostly translated the whole part
    # up, and targets left from earlier rounds then pulled vertices back down: creases. The old lift is LIFT_ANCHOR=0,
    # LIFT_ONE_SIDED=false, LIFT_HOLD=0, LIFT_PRECOND_FLOOR=-1.
    LIFT_ANCHOR=1.0,             # > 0: while lifting, hold the grounded vertices at their height (weight x LIFT_WEIGHT)
    LIFT_ONE_SIDED=True,         # lift targets only push up, never hold a vertex down
    LIFT_HOLD=0.2,               # > 0: hold every vertex that was never in a pit at its fit-only height with a
                                 # two-sided spring (weight x LIFT_WEIGHT), so a pit's rim cannot be dragged up
    LIFT_HOLD_FALLOFF=0.0,       # mm (path length through the mesh): the hold ramps from 0 at the pit to full here
    LIFT_PRECOND_FLOOR=0.3,      # PRECOND_FLOOR for the lifting solves only (< 0: same as the fit's). Lift
                                 # springs pull the whole mesh; without a floor the preconditioned step is huge at
                                 # vertices held only by micro-tets, and a dense micro-tet cluster gets torn apart
    BARRIER_WEIGHT=0.02,         # beta
    FLIP_FREE_STAGES=10,
    FLIP_FREE_STAGE_ITERATIONS=150,
    FIT_METHOD="penalty",        # "penalty": exact linear fit, then a soft fold penalty (robust on any mesh);
                                 # "barrier": fold-free barrier + staged homotopy (strict, but can stall on fine meshes)
    FOLD_PENALTY=100.0,          # penalty: strength
    FOLD_PENALTY_DET=0.5,        # penalty: tets squashed below this volume ratio (or folded) are penalised
    PENALTY_ITERATIONS=300,
    LIFT_ITERATIONS=100,         # solver iterations per lifting round (warm-started)
    BED_PIN_WEIGHT=0.0,          # > 0: hold the part's bed face flat on the bed. Off: pinning fights the tilt near the
                                 # base (the benchy got 1475 mm of unsupported extrusion at 500)
    BED_TOL=0.3,                 # mm: vertices this close to the lowest point form the bed face (always grounded)
    SLIVER_QUALITY=0.03,         # barrier: tets with mean-ratio shape quality below this get no fold barrier
    MICRO_TET_VOLUME=1e-3,       # barrier: ... nor do tets smaller than this fraction of the median tet volume
    PRECOND_FLOOR=0.0,           # minimum per-vertex stiffness in the solver's preconditioner (x median)
)


# ---------------------------------------------------------------------------------------------- kernels

@numba.njit(cache=True, inline="always")
def _det3(m):
    return (m[0, 0] * (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
            - m[0, 1] * (m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0])
            + m[0, 2] * (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0]))


@numba.njit(cache=True, inline="always")
def _adj3(m, out):
    out[0, 0] = m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1]
    out[0, 1] = m[0, 2] * m[2, 1] - m[0, 1] * m[2, 2]
    out[0, 2] = m[0, 1] * m[1, 2] - m[0, 2] * m[1, 1]
    out[1, 0] = m[1, 2] * m[2, 0] - m[1, 0] * m[2, 2]
    out[1, 1] = m[0, 0] * m[2, 2] - m[0, 2] * m[2, 0]
    out[1, 2] = m[0, 2] * m[1, 0] - m[0, 0] * m[1, 2]
    out[2, 0] = m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0]
    out[2, 1] = m[0, 1] * m[2, 0] - m[0, 0] * m[2, 1]
    out[2, 2] = m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]


@numba.njit(cache=True, inline="always")
def _J(V, c, Dm_inv, t, out):
    """out = Ds @ Dm_inv[t] for tet t with vertex ids c."""
    a0 = V[c[0]]
    for i in range(3):
        e1 = V[c[1], i] - a0[i]; e2 = V[c[2], i] - a0[i]; e3 = V[c[3], i] - a0[i]
        for j in range(3):
            out[i, j] = e1 * Dm_inv[t, 0, j] + e2 * Dm_inv[t, 1, j] + e3 * Dm_inv[t, 2, j]


@numba.njit(parallel=True, cache=True)
def _tet_terms(V, cells, Dm_inv, vol, R, beta, guard, e_out, dDs_out):
    """Per-tet energy and dE/dDs. Returns False if any guarded tet is flipped/degenerate.
    Unguarded tets (slivers) get the rotation-fit term only: no barrier, no flip check."""
    n = cells.shape[0]
    bad = np.zeros(n, np.bool_)
    for t in numba.prange(n):
        J = np.empty((3, 3)); A = np.empty((3, 3))
        _J(V, cells[t], Dm_inv, t, J)
        if not guard[t]:
            ef = 0.0
            for i in range(3):
                for j in range(3):
                    df = J[i, j] - R[t, i, j]; ef += df * df
            e_out[t] = vol[t] * ef
            for i in range(3):
                for k in range(3):
                    s = 0.0
                    for j in range(3):
                        s += vol[t] * 2.0 * (J[i, j] - R[t, i, j]) * Dm_inv[t, k, j]
                    dDs_out[t, i, k] = s
            continue
        d = _det3(J)
        if d <= 0.0:
            bad[t] = True
            continue
        _adj3(J, A)  # J^-1 = A / d
        ef = 0.0; ej = 0.0; ei = 0.0
        for i in range(3):
            for j in range(3):
                df = J[i, j] - R[t, i, j]; ef += df * df; ej += J[i, j] * J[i, j]
                inv = A[i, j] / d; ei += inv * inv
        e_out[t] = vol[t] * (ef + beta * (ej + ei - 6.0))
        # dJ = 2 (J - R) + beta (2 J - 2 J^-T J^-1 J^-T)
        # M = Ji^T Ji Ji^T with Ji = A / d
        T1 = np.empty((3, 3)); M = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                s = 0.0
                for k in range(3):
                    s += A[k, i] * A[k, j]      # (Ji^T Ji)[i, j] * d^2
                T1[i, j] = s
        for i in range(3):
            for j in range(3):
                s = 0.0
                for k in range(3):
                    s += T1[i, k] * A[j, k]     # (Ji^T Ji Ji^T)[i, j] * d^3
                M[i, j] = s / (d * d * d)
        G = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                G[i, j] = vol[t] * (2.0 * (J[i, j] - R[t, i, j]) + beta * (2.0 * J[i, j] - 2.0 * M[i, j]))
        # dE/dDs = G @ Dm_inv^T
        for i in range(3):
            for k in range(3):
                s = 0.0
                for j in range(3):
                    s += G[i, j] * Dm_inv[t, k, j]
                dDs_out[t, i, k] = s
    return not bad.any()


@numba.njit(cache=True)
def _scatter(cells, dDs, n_points):
    g = np.zeros((n_points, 3))
    for t in range(cells.shape[0]):
        for i in range(3):
            g[cells[t, 1], i] += dDs[t, i, 0]
            g[cells[t, 2], i] += dDs[t, i, 1]
            g[cells[t, 3], i] += dDs[t, i, 2]
            g[cells[t, 0], i] -= dDs[t, i, 0] + dDs[t, i, 1] + dDs[t, i, 2]
    return g


@numba.njit(parallel=True, cache=True)
def _max_step(V, D, cells, Dm_inv, guard, tcap, rho):
    """Per guarded tet: first t in (0, tcap] where det(J(V + t D)) falls to rho * its current value (tcap if never).
    rho = 0 is the flip point; rho > 0 keeps every step from squashing any tet by more than that factor."""
    n = cells.shape[0]
    out = np.empty(n)
    for t in numba.prange(n):
        if not guard[t]:
            out[t] = tcap
            continue
        A = np.empty((3, 3)); B = np.empty((3, 3)); adjA = np.empty((3, 3)); adjB = np.empty((3, 3))
        _J(V, cells[t], Dm_inv, t, A)
        _J(D, cells[t], Dm_inv, t, B)
        _adj3(A, adjA); _adj3(B, adjB)
        c0 = _det3(A); c3 = _det3(B)
        c1 = 0.0; c2 = 0.0
        for i in range(3):
            for j in range(3):
                c1 += adjA[i, j] * B[j, i]
                c2 += adjB[i, j] * A[j, i]
        # find the first sign change of p(s) = c0 + c1 s + c2 s^2 + c3 s^3 on (0, tcap], then bisect
        best = tcap
        steps = 32
        floor = rho * c0
        prev_s = 0.0
        for k in range(1, steps + 1):
            s = tcap * k / steps
            p = c0 + s * (c1 + s * (c2 + s * c3))
            if p <= floor:
                lo = prev_s; hi = s
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    pm = c0 + mid * (c1 + mid * (c2 + mid * c3))
                    if pm > floor:
                        lo = mid
                    else:
                        hi = mid
                best = lo
                break
            prev_s = s
        out[t] = best
    return out


# --------------------------------------------------------------------------------------------- solver

def mean_ratio(P, cells):
    """Tet shape quality: 1 for a regular tet, -> 0 for a sliver (12 (3V)^(2/3) / sum of squared edge lengths)."""
    a, b, c, d = (P[cells[:, k]] for k in range(4))
    V = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)) / 6
    L2 = sum(np.sum((x - y) ** 2, axis=1) for x, y in ((a, b), (a, c), (a, d), (b, c), (b, d), (c, d)))
    return 12.0 * (3.0 * V) ** (2.0 / 3.0) / np.maximum(L2, 1e-300)


class FitProblem:
    def __init__(self, P0, cells, beta, sliver_quality=0.0, step_rho=0.0, micro_volume=0.0, precond_floor=0.0):
        self.P0 = np.ascontiguousarray(P0, dtype=np.float64)
        self.cells = np.ascontiguousarray(cells, dtype=np.int64)
        a = self.P0[self.cells[:, 0]]
        Dm = np.stack([self.P0[self.cells[:, k]] - a for k in (1, 2, 3)], axis=2)
        self.Dm_inv = np.ascontiguousarray(np.linalg.inv(Dm))
        self.vol = np.abs(np.linalg.det(Dm)) / 6
        self.beta = float(beta)
        self.step_rho = float(step_rho)
        # slivers (mean ratio below sliver_quality) would pin the whole solve: a sub-micron move flips them
        # micro-tets (tiny volume, from tiny STL triangles) likewise; together they hold a negligible share of the part
        self.guard = np.ascontiguousarray((mean_ratio(self.P0, self.cells) >= sliver_quality)
                                          & (self.vol >= micro_volume * np.median(self.vol)))
        self.n = len(self.P0)
        self.e = np.empty(len(self.cells))
        self.dDs = np.empty((len(self.cells), 3, 3))
        self.R = None
        self.lift_idx = np.zeros(0, np.int64); self.lift_t = np.zeros(0); self.lift_w = np.zeros(0)
        self.lift_up = np.zeros(0, bool)
        # stiffness (Hessian of sum vol |J|^2 per coordinate, /2)
        Gi = self.Dm_inv  # rows of Dm^-1 are the basis-function gradients (J = Ds Dm^-1)
        G = np.concatenate([-Gi.sum(axis=1, keepdims=True), Gi], axis=1)
        Ke = self.vol[:, None, None] * np.einsum("nki,nli->nkl", G, G)
        c = self.cells
        K = coo_matrix((Ke.ravel(), (np.repeat(c, 4, axis=1).ravel(), np.tile(c, (1, 4)).ravel())),
                       shape=(self.n, self.n)).tocsc()
        # Vertices attached only through near-flat micro-tets have almost no stiffness, and the preconditioned step
        # blows up there (then the fold-free cap shrinks every step to ~nothing). Floor each vertex's stiffness.
        self.K_raw = K
        self.set_precond_floor(precond_floor)

    def set_precond_floor(self, precond_floor):
        dK = self.K_raw.diagonal()
        floor = precond_floor * np.median(dK)
        self.K = (self.K_raw + diags(np.maximum(floor - dK, 0.0) + 1e-8 * np.median(dK))).tocsc()

    def energy_grad(self, x):
        V = x.reshape(-1, 3)
        if not _tet_terms(V, self.cells, self.Dm_inv, self.vol, self.R, self.beta, self.guard, self.e, self.dDs):
            return np.inf, None
        E = float(np.sum(self.e))
        g = _scatter(self.cells, self.dDs, self.n)
        return E + self._lift_terms(V, g), g.ravel()

    def _lift_terms(self, V, g):
        """Soft z targets: w (z - t)^2, or only while z < t for the one-sided ones. Adds to g, returns the energy."""
        if not len(self.lift_idx):
            return 0.0
        dz = V[self.lift_idx, 2] - self.lift_t
        if self.lift_up.any():
            dz = np.where(self.lift_up & (dz > 0.0), 0.0, dz)
        np.add.at(g[:, 2], self.lift_idx, 2.0 * self.lift_w * dz)
        return float(np.sum(self.lift_w * dz * dz))

    def max_step(self, x, d):
        return float(_max_step(x.reshape(-1, 3), np.ascontiguousarray(d.reshape(-1, 3)), self.cells, self.Dm_inv, self.guard, 1.2, self.step_rho).min())

    def preconditioner(self):
        scale = 2.0 * (1.0 + 2.0 * self.beta)
        Kx = splu((scale * self.K).tocsc())
        wz = np.zeros(self.n)
        if len(self.lift_idx):
            np.add.at(wz, self.lift_idx, 2.0 * self.lift_w)
        Kz = splu((scale * self.K + diags(wz)).tocsc()) if len(self.lift_idx) else Kx

        def H0(q):
            Q = q.reshape(-1, 3)
            return np.column_stack([Kx.solve(np.ascontiguousarray(Q[:, 0])), Kx.solve(np.ascontiguousarray(Q[:, 1])),
                                    Kz.solve(np.ascontiguousarray(Q[:, 2]))]).ravel()
        return H0

    def solve(self, R, x0, max_iter, lift=None, m=10, tol=1e-7):
        self.R = np.ascontiguousarray(R)
        if lift is None:
            lift = (np.zeros(0, np.int64), np.zeros(0), np.zeros(0))
        self.lift_idx, self.lift_t, self.lift_w = lift[:3]
        # optional 4th entry: True = one-sided target (only pushes up, never holds a vertex down)
        self.lift_up = np.asarray(lift[3], bool) if len(lift) > 3 else np.zeros(len(self.lift_idx), bool)
        H0 = self.preconditioner()
        x = x0.ravel().copy()
        E, g = self.energy_grad(x)
        if not np.isfinite(E):
            raise RuntimeError("fold-free fit: starting point has inverted tets")
        S, Y = [], []
        it = 0
        for it in range(max_iter):
            q = g.copy(); al = []
            for s, y in zip(reversed(S), reversed(Y)):
                a = (s @ q) / (y @ s); al.append(a); q -= a * y
            q = H0(q)
            for (s, y), a in zip(zip(S, Y), reversed(al)):
                b = (y @ q) / (y @ s); q += s * (a - b)
            d = -q
            if g @ d >= 0:
                d = -H0(g); S, Y = [], []
            t = min(1.0, (1.0 if self.step_rho > 0 else 0.9) * self.max_step(x, d))
            for _ in range(40):
                xn = x + t * d
                En, gn = self.energy_grad(xn)
                if np.isfinite(En) and En <= E + 1e-4 * t * (g @ d):
                    break
                t *= 0.5
            else:
                break
            s = xn - x; y = gn - g
            if y @ s > 1e-12:
                S.append(s); Y.append(y)
                if len(S) > m:
                    S.pop(0); Y.pop(0)
            rel = (E - En) / max(abs(E), 1e-12)
            x, E, g = xn, En, gn
            if rel < tol and it > 5:
                break
        return x.reshape(-1, 3), it + 1


# ---------------------------------------------------------------------------------------- island tools

def vertex_graph(cells, n):
    e = np.vstack([cells[:, [i, j]] for i in range(4) for j in range(i + 1, 4)])
    e = np.unique(np.sort(e, axis=1), axis=0)
    nbrs = [[] for _ in range(n)]
    for a, b in e.tolist():
        nbrs[a].append(b); nbrs[b].append(a)
    return nbrs


def priority_flood(V, nbrs, slope, bed_tol, bed=None):
    """Lowest heights >= z so that every vertex is reachable from the bed by a path rising >= slope per mm. The bed
    is the given vertices (the pinned bed face) or else everything within bed_tol of the lowest point."""
    z = V[:, 2]; n = len(z)
    h = np.full(n, np.inf); done = np.zeros(n, bool); pq = []
    for v in (np.nonzero(z <= z.min() + bed_tol)[0] if bed is None else bed).tolist():
        h[v] = z[v]; heapq.heappush(pq, (float(z[v]), v))
    xy = V[:, :2].tolist(); zl = z.tolist()
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
            cand = max(zl[u], hv + slope * ((xu - xv) ** 2 + (yu - yv) ** 2) ** 0.5)
            if cand < h[u]:
                h[u] = cand; heapq.heappush(pq, (cand, u))
    return h


def island_seeds(V, nbrs, bed_tol=1.0, min_persistence=0.2, bed=None):
    """Local height minima not connected to the bed, with how much height they float for (sublevel-set persistence).
    Returns [(birth vertex, birth z, merge z, vertices)]."""
    z = V[:, 2]; n = len(z)
    parent = list(range(n)); birth = z.copy(); size = np.ones(n, int)
    if bed is None:
        grounded = z <= z.min() + bed_tol
    else:
        grounded = np.zeros(n, bool); grounded[bed] = True

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    active = np.zeros(n, bool); seeds = []
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
            parent[young] = old; size[old] += size[young]
            grounded[old] = grounded[old] or grounded[young]; birth[old] = min(birth[old], birth[young])
    return seeds


# ----------------------------------------------------------------------------------------------- driver

def _grounded(V, bed, tol=0.5):
    """Vertices standing on the bed: the part's bed face, plus anything within tol of the lowest point."""
    return np.union1d(bed, np.nonzero(V[:, 2] <= V[:, 2].min() + tol)[0]).astype(np.int64)


def _hold_targets(V, nbrs, pit, held, weight, falloff):
    """Two-sided z springs (per-vertex weight) at V's height for the held vertices outside the pit, ramping (smoothstep) from 0 at
    the pit to full weight at a path length of falloff mm through the mesh."""
    keep = held.copy(); keep[pit] = False
    idx = np.nonzero(keep)[0].astype(np.int64)
    w = weight[idx].copy()
    if falloff > 0 and len(pit):
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import dijkstra
        a = np.repeat(np.arange(len(nbrs)), [len(x) for x in nbrs]); b = np.concatenate([np.asarray(x, np.int64) for x in nbrs])
        G = csr_matrix((np.linalg.norm(V[a] - V[b], axis=1) + 1e-9, (a, b)), shape=(len(nbrs), len(nbrs)))
        d = dijkstra(G, directed=False, indices=pit, min_only=True, limit=falloff)
        t = np.clip(d[idx] / falloff, 0.0, 1.0)
        w *= t * t * (3.0 - 2.0 * t)
    return idx, V[idx, 2].copy(), w, np.zeros(len(idx), bool)


def deform(points, cells, cell_centers, rotation_field, p, lift=True, log=None):
    """Fold-free fit of the rotation field, then island lifting. Returns (new vertex positions, info dict)."""
    q = {**DEFAULTS, **{k: p[k] for k in DEFAULTS if k in p}}
    t0 = time.perf_counter()
    P0 = np.asarray(points, dtype=np.float64)
    V = P0.copy(); iters = 0
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
        for _ in range(int(q["ISLAND_LIFT_ROUNDS"])):
            h = priority_flood(V, nbrs, float(q["ISLAND_LIFT_SLOPE"]), bed_tol=0.5,
                               bed=_grounded(V, bed))
            need = (h - V[:, 2]) > 0.02
            if not need.any():
                break
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
    about_t = np.degrees(np.einsum("ij,ij->i", Rotation.from_matrix(U @ Wt).as_rotvec(), tangential_vectors(cell_centers)))
    want = np.degrees(rotation_field)
    sel = np.abs(want) > 20.0
    def _wmed(x, w):
        o = np.argsort(x); cw = np.cumsum(w[o]); return float(x[o][np.searchsorted(cw, cw[-1] / 2)])
    # tilt that matters for the B axis: rotation about the tangential axis, volume-weighted, where >20 deg is wanted
    achieved = abs(_wmed(about_t[sel], prob.vol[sel])) if sel.any() else 0.0
    wanted = abs(_wmed(want[sel], prob.vol[sel])) if sel.any() else 0.0
    info = {"seconds": round(time.perf_counter() - t0, 2), "fit_seconds": round(t_fit, 2), "iterations": iters,
            "inverted_tets": int((det <= 0).sum()),
            "inverted_volume_pct": round(100.0 * float(prob.vol[det <= 0].sum() / prob.vol.sum()), 4),
            "tilt_deg": round(achieved, 1), "target_tilt_deg": round(wanted, 1),
            "lifted_vertices": lifted, "island_seeds": len(seeds),
            "island_mass": round(float(sum((d - b) * m for _, b, d, m in seeds)), 1)}
    if log:
        log(f"[deform] island-free: {info}")
        if wanted > 5 and achieved < 0.7 * wanted:
            log(f"[deform] WARNING: the deformation reached only {achieved:.0f} deg of the {wanted:.0f} deg tilt it aimed "
                f"for (volume-weighted, where >20 deg is wanted). The nozzle will tilt too little; try FIT_METHOD \"penalty\".")
    return V, info


# ---------------------------------------------------------------------- penalty (soft anti-fold) solver

@numba.njit(parallel=True, cache=True)
def _tet_terms_penalty(V, cells, Dm_inv, vol, R, gamma, eps, e_out, dDs_out):
    """Per-tet energy vol * (|J - R|^2 + gamma * max(0, eps - det J)^2) and dE/dDs. Never infinite: folded tets
    are pushed back open by the penalty instead of being forbidden (so no step-size cap is needed)."""
    n = cells.shape[0]
    for t in numba.prange(n):
        J = np.empty((3, 3)); A = np.empty((3, 3))
        _J(V, cells[t], Dm_inv, t, J)
        d = _det3(J)
        _adj3(J, A)               # adj(J); d det / dJ = adj(J)^T
        viol = eps - d
        ef = 0.0
        for i in range(3):
            for j in range(3):
                df = J[i, j] - R[t, i, j]; ef += df * df
        pen = viol * viol if viol > 0.0 else 0.0
        e_out[t] = vol[t] * (ef + gamma * pen)
        G = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                g = 2.0 * (J[i, j] - R[t, i, j])
                if viol > 0.0:
                    g -= 2.0 * gamma * viol * A[j, i]
                G[i, j] = vol[t] * g
        for i in range(3):
            for k in range(3):
                s = 0.0
                for j in range(3):
                    s += G[i, j] * Dm_inv[t, k, j]
                dDs_out[t, i, k] = s


class PenaltyProblem(FitProblem):
    """Same fit as FitProblem, with a soft fold penalty instead of the barrier (robust on meshes with slivers)."""

    def __init__(self, P0, cells, gamma=100.0, eps=0.2, precond_floor=0.2):
        super().__init__(P0, cells, 0.0, precond_floor=precond_floor)
        self.gamma = float(gamma); self.eps = float(eps)

    def energy_grad(self, x):
        V = x.reshape(-1, 3)
        _tet_terms_penalty(V, self.cells, self.Dm_inv, self.vol, self.R, self.gamma, self.eps, self.e, self.dDs)
        E = float(np.sum(self.e))
        g = _scatter(self.cells, self.dDs, self.n)
        return E + self._lift_terms(V, g), g.ravel()

    def max_step(self, x, d):
        return 1.2  # no hard constraint

    def linear_start(self, R, pin=None):
        """Exact minimiser of the fit term (sum vol |J - R|^2), plus the z targets in pin if given: the starting
        point."""
        Gi = self.Dm_inv  # rows of Dm^-1 are the basis-function gradients (J = Ds Dm^-1)
        G = np.concatenate([-Gi.sum(axis=1, keepdims=True), Gi], axis=1)
        T = self.vol[:, None, None] * R
        B = np.einsum("nij,nkj->nki", T, G)
        b = np.zeros((self.n, 3))
        for k in range(4):
            np.add.at(b, self.cells[:, k], B[:, k, :])
        K = self.K_raw
        reg = 1e-9 * np.median(K.diagonal())
        lu = splu((K + reg * identity(self.n, format="csc")).tocsc())
        cols = [lu.solve(b[:, k] + 1e-9 * self.P0[:, k]) for k in range(3)]
        if pin is not None:
            idx, t, w = pin
            wz = np.zeros(self.n); bz = b[:, 2] + 1e-9 * self.P0[:, 2]
            np.add.at(wz, idx, w); np.add.at(bz, idx, w * t)
            cols[2] = splu((K + diags(wz + reg)).tocsc()).solve(bz)
        return np.column_stack(cols)
