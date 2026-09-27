"""
Island-free, fold-free deformation (DEFORMATION_METHOD = "island_free").

The notebook turns the optimised rotation field into a deformed mesh with an unconverged least-squares solve.
That solve folds the mesh (202 inverted tets on the benchy) and, more importantly, leaves local height minima
in the deformed shape. A planar slicer starts printing every such minimum in mid-air: those are the floating
islands. This module replaces only that step; the rotation field is still the notebook's.

1. Fold-free fit. Minimise over vertex positions V
       sum_c vol_c ( |J_c - R_c|^2  +  beta (|J_c|^2 + |J_c^-1|^2 - 6) )
   J_c = deformation gradient of tet c, R_c = its target rotation. The second term (symmetric Dirichlet) is a
   barrier that is infinite at zero volume, and every step is capped before the first tet would flip, so the mesh
   can never fold. Large rotations are reached by ramping the targets up in stages (homotopy). L-BFGS,
   preconditioned with the mesh stiffness matrix.
2. Island lifting. Priority-flood from the bed gives, for every vertex, the lowest height at which it can be
   reached from the bed by a path rising at least ISLAND_LIFT_SLOPE per mm. Vertices below that height sit in a
   "pit" that would print in mid-air; they get soft height targets and the fold-free fit is re-solved, a few rounds.

Deterministic: all per-tet work is written to arrays and reduced serially.
"""
import heapq
import time

import numba
import numpy as np
from scipy.sparse import coo_matrix, diags, identity
from scipy.sparse.linalg import splu

from .fast_map import rotation_matrices

DEFAULTS = dict(
    ISLAND_LIFT_SLOPE=0.5,       # mm of rise per mm; 0.5 ~ overhangs up to ~63 deg from vertical in deformed space
    ISLAND_LIFT_ROUNDS=5,
    LIFT_WEIGHT=5.0,
    BARRIER_WEIGHT=0.02,         # beta
    FLIP_FREE_STAGES=10,
    FLIP_FREE_STAGE_ITERATIONS=150,
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
def _tet_terms(V, cells, Dm_inv, vol, R, beta, e_out, dDs_out):
    """Per-tet energy and dE/dDs. Returns False if any tet is flipped/degenerate."""
    n = cells.shape[0]
    bad = np.zeros(n, np.bool_)
    for t in numba.prange(n):
        J = np.empty((3, 3)); A = np.empty((3, 3))
        _J(V, cells[t], Dm_inv, t, J)
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
def _max_step(V, D, cells, Dm_inv, tcap):
    """Per tet: first t in (0, tcap] where det(J(V + t D)) hits 0 (tcap if none)."""
    n = cells.shape[0]
    out = np.empty(n)
    for t in numba.prange(n):
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
        prev_s = 0.0; prev_p = c0
        for k in range(1, steps + 1):
            s = tcap * k / steps
            p = c0 + s * (c1 + s * (c2 + s * c3))
            if p <= 0.0:
                lo = prev_s; hi = s
                for _ in range(40):
                    mid = 0.5 * (lo + hi)
                    pm = c0 + mid * (c1 + mid * (c2 + mid * c3))
                    if pm > 0.0:
                        lo = mid
                    else:
                        hi = mid
                best = lo
                break
            prev_s = s; prev_p = p
        out[t] = best
    return out


# --------------------------------------------------------------------------------------------- solver

class FitProblem:
    def __init__(self, P0, cells, beta):
        self.P0 = np.ascontiguousarray(P0, dtype=np.float64)
        self.cells = np.ascontiguousarray(cells, dtype=np.int64)
        a = self.P0[self.cells[:, 0]]
        Dm = np.stack([self.P0[self.cells[:, k]] - a for k in (1, 2, 3)], axis=2)
        self.Dm_inv = np.ascontiguousarray(np.linalg.inv(Dm))
        self.vol = np.abs(np.linalg.det(Dm)) / 6
        self.beta = float(beta)
        self.n = len(self.P0)
        self.e = np.empty(len(self.cells))
        self.dDs = np.empty((len(self.cells), 3, 3))
        self.R = None
        self.lift_idx = np.zeros(0, np.int64); self.lift_t = np.zeros(0); self.lift_w = np.zeros(0)
        # stiffness (Hessian of sum vol |J|^2 per coordinate, /2)
        Gi = np.swapaxes(self.Dm_inv, 1, 2)
        G = np.concatenate([-Gi.sum(axis=1, keepdims=True), Gi], axis=1)
        Ke = self.vol[:, None, None] * np.einsum("nki,nli->nkl", G, G)
        c = self.cells
        self.K = coo_matrix((Ke.ravel(), (np.repeat(c, 4, axis=1).ravel(), np.tile(c, (1, 4)).ravel())),
                            shape=(self.n, self.n)).tocsc() + 1e-8 * identity(self.n, format="csc")

    def energy_grad(self, x):
        V = x.reshape(-1, 3)
        if not _tet_terms(V, self.cells, self.Dm_inv, self.vol, self.R, self.beta, self.e, self.dDs):
            return np.inf, None
        E = float(np.sum(self.e))
        g = _scatter(self.cells, self.dDs, self.n)
        if len(self.lift_idx):
            dz = V[self.lift_idx, 2] - self.lift_t
            E += float(np.sum(self.lift_w * dz * dz))
            np.add.at(g[:, 2], self.lift_idx, 2.0 * self.lift_w * dz)
        return E, g.ravel()

    def max_step(self, x, d):
        return float(_max_step(x.reshape(-1, 3), np.ascontiguousarray(d.reshape(-1, 3)), self.cells, self.Dm_inv, 1.2).min())

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
            self.lift_idx = np.zeros(0, np.int64); self.lift_t = np.zeros(0); self.lift_w = np.zeros(0)
        else:
            self.lift_idx, self.lift_t, self.lift_w = lift
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
            t = min(1.0, 0.9 * self.max_step(x, d))
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


def priority_flood(V, nbrs, slope, bed_tol):
    """Lowest heights >= z so that every vertex is reachable from the bed by a path rising >= slope per mm."""
    z = V[:, 2]; n = len(z)
    h = np.full(n, np.inf); done = np.zeros(n, bool); pq = []
    for v in np.nonzero(z <= z.min() + bed_tol)[0].tolist():
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


def island_seeds(V, nbrs, bed_tol=1.0, min_persistence=0.2):
    """Local height minima not connected to the bed, with how much height they float for (sublevel-set persistence).
    Returns [(birth vertex, birth z, merge z, vertices)]."""
    z = V[:, 2]; n = len(z)
    parent = list(range(n)); birth = z.copy(); grounded = z <= z.min() + bed_tol; size = np.ones(n, int)

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

def deform(points, cells, cell_centers, rotation_field, p, lift=True, log=None):
    """Fold-free fit of the rotation field, then island lifting. Returns (new vertex positions, info dict)."""
    q = {**DEFAULTS, **{k: p[k] for k in DEFAULTS if k in p}}
    t0 = time.perf_counter()
    P0 = np.asarray(points, dtype=np.float64)
    prob = FitProblem(P0, cells, q["BARRIER_WEIGHT"])
    V = P0.copy(); iters = 0
    stages = int(q["FLIP_FREE_STAGES"])
    R = None
    for k in range(1, stages + 1):
        R = rotation_matrices(cell_centers, (k / stages) * rotation_field)
        V, it = prob.solve(R, V, int(q["FLIP_FREE_STAGE_ITERATIONS"]))
        iters += it
    t_fit = time.perf_counter() - t0
    nbrs = vertex_graph(np.asarray(cells), len(P0))
    lifted = 0
    if lift:
        tgt = np.full(len(P0), -np.inf)
        for r in range(int(q["ISLAND_LIFT_ROUNDS"])):
            h = priority_flood(V, nbrs, float(q["ISLAND_LIFT_SLOPE"]), bed_tol=0.5)
            need = (h - V[:, 2]) > 0.02
            if not need.any():
                break
            tgt = np.maximum(tgt, np.where(need, h, -np.inf))
            idx = np.nonzero(np.isfinite(tgt))[0].astype(np.int64)
            V, it = prob.solve(R, V, 300, lift=(idx, tgt[idx], np.full(len(idx), float(q["LIFT_WEIGHT"]))))
            iters += it
            lifted = len(idx)
    seeds = island_seeds(V, nbrs)
    a = V[np.asarray(cells)[:, 0]]
    det = np.linalg.det(np.stack([V[np.asarray(cells)[:, k]] - a for k in (1, 2, 3)], axis=2))
    info = {"seconds": round(time.perf_counter() - t0, 2), "fit_seconds": round(t_fit, 2), "iterations": iters,
            "inverted_tets": int((det <= 0).sum()), "lifted_vertices": lifted, "island_seeds": len(seeds),
            "island_mass": round(float(sum((d - b) * m for _, b, d, m in seeds)), 1)}
    if log:
        log(f"[deform] island-free: {info}")
    return V, info
