"""
Bit-identical acceleration of scipy's TRF least_squares for our sparse Jacobians.

scipy's TRF (tr_solver='lsmr') spends most of its time in single-threaded sparse
mat-vecs J@x and J.T@v. Those are plain per-row sums in index order starting from 0.0,
so a row-parallel numba kernel over CSR (for J) and over CSR-of-J.T (for J.T) produces
exactly the same float64 results, ~20x faster.

`fast_trf(jacobian_builder)` temporarily swaps the two helpers scipy's trf module uses to
wrap J (compute_grad, right_multiplied_operator) for versions that call the numba kernels
when J is the matrix our builder produced. Everything else in TRF (scaling via
J.power(2).sum(), LSMR itself, step control) is untouched.
"""
from contextlib import contextmanager

import numba
import numpy as np
from scipy.optimize._lsq import trf as _trf
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.linalg import LinearOperator


@numba.njit(parallel=True, fastmath=False, cache=True)
def csr_matvec(indptr, indices, data, x, out):
    for i in numba.prange(len(indptr) - 1):
        s = 0.0
        for k in range(indptr[i], indptr[i + 1]):
            s += np.float64(data[k]) * x[indices[k]]
        out[i] = s
    return out


class FastJacobian:
    """Builds CSR Jacobians with a fixed pattern (as lil_matrix+tocsr would) and remembers
    the transpose layout so J.T@v can run row-parallel as well."""

    def __init__(self, n_rows, n_cols, rows, cols):
        self.shape = (n_rows, n_cols)
        self.rows, self.cols = rows, cols
        self.order = np.lexsort((cols, rows))
        self.indices = cols[self.order].astype(np.int32)
        self.indptr = np.concatenate([[0], np.cumsum(np.bincount(rows, minlength=n_rows))]).astype(np.int32)
        # transpose: entries sorted by (col, row) -> CSR of J.T with rows ascending inside each column
        sorted_rows = rows[self.order]
        self.t_perm = np.lexsort((sorted_rows, self.indices))
        self.t_indices = sorted_rows[self.t_perm].astype(np.int32)
        self.t_indptr = np.concatenate([[0], np.cumsum(np.bincount(self.indices, minlength=n_cols))]).astype(np.int32)
        self.sorted_rows = sorted_rows                 # row of each entry in CSR order
        self.t_rows = self.indices[self.t_perm]        # row of J.T (= column of J) of each transposed entry
        self.current = None

    def build(self, vals):
        data = vals.astype(np.float32)[self.order]
        if data.all():
            J = csr_matrix((data, self.indices, self.indptr), shape=self.shape, copy=False)
            T = (self.t_indptr, self.t_indices, data[self.t_perm])
        else:  # lil_matrix never stores explicit zeros: drop them from the fixed layout
            nz = data != 0
            indptr = np.concatenate([[0], np.cumsum(np.bincount(self.sorted_rows[nz], minlength=self.shape[0]))]).astype(np.int32)
            J = csr_matrix((data[nz], self.indices[nz], indptr), shape=self.shape, copy=False)
            t_data = data[self.t_perm]
            t_nz = t_data != 0
            t_indptr = np.concatenate([[0], np.cumsum(np.bincount(self.t_rows[t_nz], minlength=self.shape[1]))]).astype(np.int32)
            T = (t_indptr, self.t_indices[t_nz], t_data[t_nz])
        self.current = (J.data, J.indptr, J.indices, T)
        return J

    def owns(self, J):
        """True if J is (a copy of) the matrix built last. scipy's VectorFunction copies the
        Jacobian into a csr_array, so compare content rather than identity."""
        if self.current is None or not issparse(J) or J.format != "csr" or J.shape != self.shape:
            return False
        data, indptr, indices, _ = self.current
        return (J.nnz == len(data) and np.array_equal(J.indptr, indptr)
                and np.array_equal(J.indices, indices) and np.array_equal(J.data, data))

    def matvec(self, x):
        data, indptr, indices, _ = self.current
        return csr_matvec(indptr, indices, data, np.ascontiguousarray(x, dtype=np.float64), np.empty(self.shape[0]))

    def rmatvec(self, v):
        t_indptr, t_indices, t_data = self.current[3]
        return csr_matvec(t_indptr, t_indices, t_data, np.ascontiguousarray(v, dtype=np.float64), np.empty(self.shape[1]))


@contextmanager
def fast_trf(fj):
    orig_grad = _trf.compute_grad
    orig_rmo = _trf.right_multiplied_operator

    def compute_grad(J, f):
        if fj.owns(J):
            return fj.rmatvec(f)  # == J.T.dot(f)
        return orig_grad(J, f)

    def right_multiplied_operator(J, d):
        if not fj.owns(J):
            return orig_rmo(J, d)

        def matvec(x):
            return fj.matvec(np.ravel(x) * d)

        def matmat(X):
            Xd = X * d[:, np.newaxis]
            return np.column_stack([fj.matvec(Xd[:, k]) for k in range(Xd.shape[1])])

        def rmatvec(x):
            return d * fj.rmatvec(np.ravel(x))

        return LinearOperator(J.shape, matvec=matvec, matmat=matmat, rmatvec=rmatvec, dtype=np.float64)

    _trf.compute_grad = compute_grad
    _trf.right_multiplied_operator = right_multiplied_operator
    try:
        yield
    finally:
        _trf.compute_grad = orig_grad
        _trf.right_multiplied_operator = orig_rmo
