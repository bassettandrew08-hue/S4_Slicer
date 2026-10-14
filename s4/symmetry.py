"""
Mirror symmetry for (nearly) mirror-symmetric parts (SYMMETRY_PLANE, not in the notebook; off by default).

TetGen's mesh of a symmetric part is not symmetric, so neither is the rotation (tilt) field optimised on it nor the
deformed shape, and the two halves of the part print with different curved layers (pi 3mm: the fit leans the whole
part ~1 deg about y, 0.74 mm of layer mismatch between its legs). With a mirror plane (vertical, through the part's
surface centroid; SYMMETRY_PLANE = its normal's angle in the xy plane in deg, a list of angles, or "auto" to search
for them; with two planes the averages are applied one after the other) this module averages, with weight w = 1/2:

  "field"    each tet's rotation with the rotation at its mirrored centroid. The tilt is about the tangential axis,
             and a mirror maps a rotation about the tangential axis at p to the same angle about the tangential axis
             at the mirrored point (for a plane through the bed's z axis), so the scalar angle is averaged as it is.
             This changes Joshua Bird's field (only when the setting is on).
  "vertices" each deformed vertex with the mirror image of the deformation at its mirrored point (barycentric
             interpolation in the undeformed mesh), mirrored across the deformed shape's own plane (same normal,
             offset fitted). Applied between island_free's fit and its lift, so the lift still removes any pit the
             average creates. Weights are halved where the average would squash or fold a tet. Any fold-free
             deformation prints the original part, so this only changes the layer shapes.

Where the mirrored point lies outside the part (details that are not mirrored, such as text) w ramps from 1/2 at
the surface to 0 at SYMMETRY_TOLERANCE mm outside, so those regions keep their own values.
"""
import numpy as np

# "auto" accepts the best plane only if the mirrored surface lies within this median distance (mm) of the surface:
# dino 0.08, benchy 0.07, pi 0.0; Squirtle, z mount and B axis mount (not symmetric) 1.1
AUTO_MAX_ERROR = 0.25


def surface_samples(points, cells, n=60000, seed=0):
    """(triangle mesh of the tet mesh's boundary, uniform samples on it)."""
    import open3d as o3d
    f = np.vstack([cells[:, [0, 1, 2]], cells[:, [0, 1, 3]], cells[:, [0, 2, 3]], cells[:, [1, 2, 3]]])
    key = np.sort(f, axis=1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    tri = f[cnt[inv.reshape(-1)] == 1]
    m = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(np.asarray(points, float)),
                                  o3d.utility.Vector3iVector(tri))
    m.remove_unreferenced_vertices()
    np.random.seed(seed)
    o3d.utility.random.seed(seed)
    s = np.asarray(m.sample_points_uniformly(n).points)
    return m, s


def _scene(mesh):
    import open3d as o3d
    sc = o3d.t.geometry.RaycastingScene()
    sc.add_triangles(o3d.t.geometry.TriangleMesh.from_legacy(mesh))
    return sc


def _dist(scene, P, signed=False):
    import open3d as o3d
    T = o3d.core.Tensor(np.asarray(P, dtype=np.float32))
    return (scene.compute_signed_distance(T) if signed else scene.compute_distance(T)).numpy().astype(float)


def mirror(P, n, d):
    """Reflect points across the plane x.n = d."""
    return P - 2.0 * ((P @ n) - d)[:, None] * n


def normal(angle_deg):
    a = np.radians(float(angle_deg))
    return np.array([np.cos(a), np.sin(a), 0.0])


def find_planes(points, cells, setting, log=print):
    """-> [(normal, offset, info)] of the mirror planes (vertical, through the surface centroid). setting: an angle
    in deg, a list of angles, or "auto" (every plane whose mirrored surface is within AUTO_MAX_ERROR, median, of the
    surface; pi 3mm has two)."""
    mesh, S = surface_samples(points, cells)
    V, F = np.asarray(mesh.vertices), np.asarray(mesh.triangles)
    area = np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1)
    c = (V[F].mean(1) * area[:, None]).sum(0) / area.sum()  # exact surface centroid (no sampling noise)
    sc = _scene(mesh)
    sub = S[::6]

    def score(a):
        n = normal(a)
        return float(np.median(_dist(sc, mirror(sub, n, c @ n))))
    if isinstance(setting, str):
        if setting.lower() != "auto":
            raise ValueError(f"SYMMETRY_PLANE must be null, an angle in deg, a list of angles or \"auto\", "
                             f"not {setting!r}")
        angles = np.arange(0.0, 180.0, 2.0)
        s = np.array([score(a) for a in angles])
        angles_found = []
        for i in np.argsort(s):
            if s[i] > 2 * AUTO_MAX_ERROR or s[i] > min(s[i - 1], s[(i + 1) % len(s)]):
                continue  # not a local minimum (circular in the angle) or hopeless
            best = angles[i]
            for step in (0.5, 0.1):
                fine = best + step * np.arange(-4, 5)
                best = fine[np.argmin([score(a) for a in fine])]
            best = float(best % 180.0)
            if all(min(abs(best - a), 180 - abs(best - a)) > 5 for a in angles_found):
                angles_found.append(best)
        if not angles_found:
            angles_found = [float(angles[np.argmin(s)])]
    else:
        angles_found = [float(a) for a in (setting if isinstance(setting, (list, tuple)) else [setting])]
    planes = []
    for angle in angles_found:
        n = normal(angle)
        err = _dist(sc, mirror(S, n, c @ n))
        info = {"angle_deg": round(angle, 2), "offset_mm": round(float(c @ n), 3),
                "mirror_err_median_mm": round(float(np.median(err)), 3),
                "mirror_err_p90_mm": round(float(np.percentile(err, 90)), 3)}
        ok = not isinstance(setting, str) or np.median(err) <= AUTO_MAX_ERROR
        log(f"[symmetry] mirror plane: normal {angle:.1f} deg in xy, through the surface centroid "
            f"(x.n = {c @ n:.3f}); mirrored surface off the surface by median {info['mirror_err_median_mm']} mm, "
            f"90th pct {info['mirror_err_p90_mm']} mm" + ("" if ok else f": over {AUTO_MAX_ERROR} mm, not used"))
        if ok:
            planes.append((n, float(c @ n), info))
    if not planes:
        log("[symmetry] no mirror plane found; symmetry off")
    return planes


def _locate(grid, cells, points, Q):
    """Containing (else closest) tet of each query point and its barycentric coordinates (clipped to >= 0)."""
    ci = np.asarray(grid.find_containing_cell(Q))
    miss = ci < 0
    if miss.any():
        ci[miss] = np.asarray(grid.find_closest_cell(Q[miss]))
    T = points[cells[ci]]
    A = (T[:, 1:] - T[:, :1]).transpose(0, 2, 1)
    lam = np.linalg.solve(A, (Q - T[:, 0])[:, :, None])[:, :, 0]
    b = np.c_[1.0 - lam.sum(1), lam]
    b = np.clip(b, 0.0, None)
    b /= b.sum(1, keepdims=True)
    return ci, b


class _Plane:
    pass


class Mirror:
    """The planes plus the mirrored-point lookups on the undeformed mesh (computed once, reused per iteration).
    With two planes the averages are applied one after the other."""

    def __init__(self, grid, cells, setting, tolerance, log=print):
        self.grid, self.cells = grid, np.asarray(cells)
        self.points = np.asarray(grid.points, dtype=float)
        self.log = log
        self.tol = float(tolerance)
        self.planes = []
        found = find_planes(self.points, self.cells, setting, log=log)
        self.ok = bool(found)
        if not self.ok:
            return
        mesh, _ = surface_samples(self.points, self.cells, n=10)
        sc = _scene(mesh)

        def weight(Q):  # 1/2 inside, ramping to 0 at tol outside
            sd = _dist(sc, Q, signed=True)
            out = np.clip(sd, 0.0, None)
            return 0.5 * np.clip(1.0 - out / max(self.tol, 1e-9), 0.0, 1.0)
        cc = np.asarray(grid.cell_centers().points)
        for n, d, info in found:
            pl = _Plane()
            pl.n, pl.d, pl.info = n, d, info
            cm = mirror(cc, n, d)
            pl.cell_w = weight(cm)
            pl.cell_ci, pl.cell_b = _locate(grid, self.cells, self.points, cm)
            vm = mirror(self.points, n, d)
            pl.vert_w = weight(vm)
            pl.vert_ci, pl.vert_b = _locate(grid, self.cells, self.points, vm)
            self.planes.append(pl)

    def field(self, rf):
        for pl in self.planes:
            rf = self._field(pl, rf)
        return rf

    def vertices(self, V):
        for pl in self.planes:
            V = self._vertices(pl, V)
        return V

    def _field(self, pl, rf):
        """Average each tet's rotation with the rotation at its mirrored centroid (piecewise-linear: tet values
        averaged to the vertices, volume-weighted, then interpolated)."""
        vol = np.abs(np.asarray(self.grid.compute_cell_sizes(length=False, area=False)["Volume"]))
        nv = len(self.points)
        num = np.zeros(nv)
        den = np.zeros(nv)
        for k in range(4):
            np.add.at(num, self.cells[:, k], vol * rf)
            np.add.at(den, self.cells[:, k], vol)
        rv = num / np.maximum(den, 1e-30)
        at_m = np.einsum("ij,ij->i", pl.cell_b, rv[self.cells[pl.cell_ci]])
        out = (1.0 - pl.cell_w) * rf + pl.cell_w * at_m
        self.log(f"[symmetry] field: {100 * np.mean(pl.cell_w > 0):.0f}% of tets averaged with their mirror; mean "
                 f"change {np.degrees(np.mean(np.abs(out - rf))):.2f} deg, max {np.degrees(np.abs(out - rf).max()):.2f}")
        return out

    def _vertices(self, pl, V):
        """Average each deformed vertex with the mirrored deformation at its mirrored point. Where that would fold or
        squash a tet (volume ratio below the input's own minimum or 0.1), the weight of its vertices is halved,
        repeatedly, so the result has no more folded tets than the input."""
        Dm = np.einsum("ij,ijk->ik", pl.vert_b, V[self.cells[pl.vert_ci]])
        w = pl.vert_w.copy()
        sel = w > 0
        d2 = float(np.median((((V + Dm) / 2) @ pl.n)[sel])) if sel.any() else float(V.mean(0) @ pl.n)
        target = mirror(Dm, pl.n, d2)
        c = self.cells

        def ratio(X):
            J = np.stack([X[c[:, k]] - X[c[:, 0]] for k in (1, 2, 3)], axis=2)
            J0 = np.stack([self.points[c[:, k]] - self.points[c[:, 0]] for k in (1, 2, 3)], axis=2)
            return np.linalg.det(J) / np.linalg.det(J0)
        floor = np.minimum(ratio(V), 0.1)
        backed = 0
        for _ in range(30):
            out = (1.0 - w)[:, None] * V + w[:, None] * target
            bad = ratio(out) < floor
            if not bad.any():
                break
            vb = np.unique(c[bad])
            backed += len(vb)
            w[vb] *= 0.5
        else:
            w[np.unique(c[ratio(out) < floor])] = 0.0
            out = (1.0 - w)[:, None] * V + w[:, None] * target
        mv = np.linalg.norm(out - V, axis=1)
        self.log(f"[symmetry] vertices: {100 * np.mean(sel):.0f}% averaged with their mirror; moved mean "
                 f"{mv.mean():.3f} mm, max {mv.max():.3f} mm ({backed} weight halvings to keep tets unfolded)")
        return out
