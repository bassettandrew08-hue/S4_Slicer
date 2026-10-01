"""
Solvers for the island-free deformation (s4/island_free.py, whose docstring explains the fit and the lift).

FitProblem: fold-free fit with a symmetric-Dirichlet barrier (FIT_METHOD "barrier"); PenaltyProblem: the same fit with
a soft fold penalty instead (FIT_METHOD "penalty", the default). Both minimise with preconditioned L-BFGS and take
soft z targets for the lift. The per-tet terms are numba kernels.

Deterministic: all per-tet work is written to arrays and reduced serially.
"""
import numba
import numpy as np
from scipy.sparse import coo_matrix, diags, identity
from scipy.sparse.linalg import splu


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
        e1 = V[c[1], i] - a0[i]
        e2 = V[c[2], i] - a0[i]
        e3 = V[c[3], i] - a0[i]
        for j in range(3):
            out[i, j] = e1 * Dm_inv[t, 0, j] + e2 * Dm_inv[t, 1, j] + e3 * Dm_inv[t, 2, j]


@numba.njit(parallel=True, cache=True)
def _tet_terms(V, cells, Dm_inv, vol, R, beta, guard, e_out, dDs_out):
    """Per-tet energy and dE/dDs. Returns False if any guarded tet is flipped/degenerate.
    Unguarded tets (slivers) get the rotation-fit term only: no barrier, no flip check."""
    n = cells.shape[0]
    bad = np.zeros(n, np.bool_)
    for t in numba.prange(n):
        J = np.empty((3, 3))
        A = np.empty((3, 3))
        _J(V, cells[t], Dm_inv, t, J)
        if not guard[t]:
            ef = 0.0
            for i in range(3):
                for j in range(3):
                    df = J[i, j] - R[t, i, j]
                    ef += df * df
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
        ef = 0.0
        ej = 0.0
        ei = 0.0
        for i in range(3):
            for j in range(3):
                df = J[i, j] - R[t, i, j]
                ef += df * df
                ej += J[i, j] * J[i, j]
                inv = A[i, j] / d
                ei += inv * inv
        e_out[t] = vol[t] * (ef + beta * (ej + ei - 6.0))
        # dJ = 2 (J - R) + beta (2 J - 2 J^-T J^-1 J^-T)
        # M = Ji^T Ji Ji^T with Ji = A / d
        T1 = np.empty((3, 3))
        M = np.empty((3, 3))
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
        A = np.empty((3, 3))
        B = np.empty((3, 3))
        adjA = np.empty((3, 3))
        adjB = np.empty((3, 3))
        _J(V, cells[t], Dm_inv, t, A)
        _J(D, cells[t], Dm_inv, t, B)
        _adj3(A, adjA)
        _adj3(B, adjB)
        c0 = _det3(A)
        c3 = _det3(B)
        c1 = 0.0
        c2 = 0.0
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
                lo = prev_s
                hi = s
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
    """Fold-free fit of a tet mesh to per-tet target rotations (FIT_METHOD 'barrier'), plus soft z targets.

    Energy over vertex positions V: sum_c vol_c (|J_c - R_c|^2 + beta (|J_c|^2 + |J_c^-1|^2 - 6)) + lift terms,
    J_c the deformation gradient of tet c (Ds Dm^-1). The barrier term is infinite at zero volume, and max_step caps
    every step before the first guarded tet would flip. Slivers and micro-tets (sliver_quality, micro_volume) are
    unguarded: rotation term only. The preconditioner is the mesh stiffness matrix K (Laplacian-like, per
    coordinate), floored at precond_floor x median per vertex. PenaltyProblem swaps the barrier for a soft penalty.
    """
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
        self.lift_idx = np.zeros(0, np.int64)
        self.lift_t = np.zeros(0)
        self.lift_w = np.zeros(0)
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
        steps = _max_step(x.reshape(-1, 3), np.ascontiguousarray(d.reshape(-1, 3)), self.cells, self.Dm_inv,
                          self.guard, 1.2, self.step_rho)
        return float(steps.min())

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
        """Minimise the energy for target rotations R (n_cells, 3, 3) from x0 with preconditioned L-BFGS (memory m).

        lift: soft z targets (vertex ids, target z, weights[, one-sided flags]); one-sided ones only push up.
        The step starts at the largest fold-free step (max_step) and is halved until the Armijo condition holds.
        Stops after max_iter iterations, when the relative energy drop falls below tol, or when no step is found.
        Returns (vertex positions (n, 3), iterations used).
        """
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
            q = g.copy()
            al = []
            for s, y in zip(reversed(S), reversed(Y)):
                a = (s @ q) / (y @ s)
                al.append(a)
                q -= a * y
            q = H0(q)
            for (s, y), a in zip(zip(S, Y), reversed(al)):
                b = (y @ q) / (y @ s)
                q += s * (a - b)
            d = -q
            if g @ d >= 0:
                d = -H0(g)
                S, Y = [], []
            t = min(1.0, (1.0 if self.step_rho > 0 else 0.9) * self.max_step(x, d))
            for _ in range(40):
                xn = x + t * d
                En, gn = self.energy_grad(xn)
                if np.isfinite(En) and En <= E + 1e-4 * t * (g @ d):
                    break
                t *= 0.5
            else:
                break
            s = xn - x
            y = gn - g
            if y @ s > 1e-12:
                S.append(s)
                Y.append(y)
                if len(S) > m:
                    S.pop(0)
                    Y.pop(0)
            rel = (E - En) / max(abs(E), 1e-12)
            x, E, g = xn, En, gn
            if rel < tol and it > 5:
                break
        return x.reshape(-1, 3), it + 1


# ---------------------------------------------------------------------- penalty (soft anti-fold) solver

@numba.njit(parallel=True, cache=True)
def _tet_terms_penalty(V, cells, Dm_inv, vol, R, gamma, eps, e_out, dDs_out):
    """Per-tet energy vol * (|J - R|^2 + gamma * max(0, eps - det J)^2) and dE/dDs. Never infinite: folded tets
    are pushed back open by the penalty instead of being forbidden (so no step-size cap is needed)."""
    n = cells.shape[0]
    for t in numba.prange(n):
        J = np.empty((3, 3))
        A = np.empty((3, 3))
        _J(V, cells[t], Dm_inv, t, J)
        d = _det3(J)
        _adj3(J, A)               # adj(J); d det / dJ = adj(J)^T
        viol = eps - d
        ef = 0.0
        for i in range(3):
            for j in range(3):
                df = J[i, j] - R[t, i, j]
                ef += df * df
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
        self.gamma = float(gamma)
        self.eps = float(eps)

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
            wz = np.zeros(self.n)
            bz = b[:, 2] + 1e-9 * self.P0[:, 2]
            np.add.at(wz, idx, w)
            np.add.at(bz, idx, w * t)
            cols[2] = splu((K + diags(wz + reg)).tocsc()).solve(bz)
        return np.column_stack(cols)
