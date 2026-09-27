"""
Fast deformation (notebook cells 2-11), numerically equivalent to s4/reference.py.

Output-preserving changes vs. the notebook:
  * neighbours: computed from connectivity with numpy, then fed through Python sets in the
    same insertion order pyvista/VTK uses, so every neighbour list has the identical order
  * networkx graph: one add_weighted_edges_from call (same edge order)
  * calculate_tet_attributes: the unused shared_vertices loop is dropped; the cell<->face
    maps use hashing instead of an O(points x surface points) search
  * multi-source Dijkstra runs once and is reused (the notebook runs the identical search 3x)
  * the in-air test uses the max height along each Dijkstra path (same comparison)
  * initial-rotation-field smoothing runs its (identical) pass once instead of 30 times:
    every pass in the notebook recomputes from the same unsmoothed field, so the result of
    pass 30 equals pass 1
  * least-squares Jacobians are assembled directly as CSR with the same float32 values and
    sparsity pattern the notebook's lil_matrix produced (zeros removed like lil does)
  * update_tet_attributes(deformed) and calculate_tet_attributes(deformed) are skipped: the
    mapper only needs connectivity, vertex positions and cell centres

The least-squares objectives and solver settings are unchanged (including the squared
residuals, i.e. the quartic penalty). See NOTES in the README before changing them.
"""
import base64
import pickle

import networkx as nx
import numpy as np
import pyvista as pv
from scipy.optimize import least_squares
from scipy.sparse import csr_matrix

from .timing import TIMER
from . import meshio_s4
from .fast_lsq import FastJacobian, fast_trf

up_vector = np.array([0, 0, 1])
N = np.eye(4) - 1 / 4 * np.ones((4, 4))

# VTK tetra local edge / face vertex orders (vtkTetra::Edges / vtkTetra::Faces)
VTK_TET_EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
VTK_TET_FACES = [(0, 1, 3), (1, 2, 3), (2, 0, 3), (0, 2, 1)]


class MeshContext:
    """Everything the notebook kept in globals for one input mesh."""
    pass


# ----------------------------------------------------------------------------- neighbours

def _group_cells_by_key(keys_per_cell):
    """keys_per_cell: (n_cells, k) int64 keys. Returns (inverse ids (n_cells,k), groups list of ascending cell lists)."""
    n_cells, k = keys_per_cell.shape
    flat = keys_per_cell.ravel()
    uniq, inv = np.unique(flat, return_inverse=True)
    inv = inv.reshape(n_cells, k)
    order = np.argsort(inv.ravel(), kind="stable")  # stable: cell-major -> ascending cells per key
    cell_of = (order // k)
    counts = np.bincount(inv.ravel(), minlength=len(uniq))
    splits = np.cumsum(counts)[:-1]
    groups = [g.tolist() for g in np.split(cell_of, splits)]
    return inv, groups


def compute_neighbours(cells, n_points):
    """
    Reproduce pyvista cell_neighbors(c, 'points'|'edges'|'faces') for every cell, including
    list order: pyvista builds a Python set by .update()-ing VTK's GetCellNeighbors results
    (ascending cell ids) for each point/edge/face in VTK's local order, then returns list(set).
    """
    n_cells = cells.shape[0]
    cells64 = cells.astype(np.int64)

    # point -> cells (ascending)
    point_inv, point_groups = _group_cells_by_key(cells64)
    # edges keyed by sorted vertex pair
    e = np.stack([np.sort(cells64[:, list(p)], axis=1) for p in VTK_TET_EDGES], axis=1)  # (n,6,2)
    edge_keys = e[:, :, 0] * n_points + e[:, :, 1]
    edge_inv, edge_groups = _group_cells_by_key(edge_keys)
    f = np.stack([np.sort(cells64[:, list(t)], axis=1) for t in VTK_TET_FACES], axis=1)  # (n,4,3)
    face_keys = (f[:, :, 0] * n_points + f[:, :, 1]) * n_points + f[:, :, 2]
    face_inv, face_groups = _group_cells_by_key(face_keys)

    result = {}
    for kind, inv, groups in (("point", point_inv, point_groups), ("edge", edge_inv, edge_groups), ("face", face_inv, face_groups)):
        inv_l = inv.tolist()
        pairs = []
        for c in range(n_cells):
            s = set()
            for gid in inv_l[c]:
                s.update([x for x in groups[gid] if x != c])
            for nb in list(s):
                if nb > c:
                    pairs.append((c, nb))
        result[kind] = np.array(pairs)
    return result


def neighbour_dicts(pairs_by_kind, n_cells, kinds=("point", "edge")):
    out = {}
    for kind in kinds:
        d = {c: [] for c in range(n_cells)}
        for a, b in pairs_by_kind[kind].tolist():
            d[a].append(b)
            d[b].append(a)
        out[kind] = d
    return out


# ------------------------------------------------------------------------ tet attributes

def cell_face_maps(tet, cells, surface_mesh):
    """Same cell_to_face / has_face as reference.calculate_tet_attributes, via hashing."""
    face_vertices = np.asarray(surface_mesh.points)
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]
    # cell vertex -> face vertex, only where the coordinate matches exactly one face vertex
    first = {}
    dup = set()
    for fi, key in enumerate(map(tuple, face_vertices.tolist())):
        if key in first:
            dup.add(key)
        else:
            first[key] = fi
    c2f_vertex = {}
    for vi, key in enumerate(map(tuple, np.asarray(tet.points).tolist())):
        if key in first and key not in dup:
            c2f_vertex[vi] = first[key]
    # surface triangles keyed by their vertex set
    tri_by_set = {}
    for ti, tri in enumerate(faces.tolist()):
        tri_by_set.setdefault(frozenset(tri), []).append(ti)
    cell_to_face = {}
    for ci, cell in enumerate(cells.tolist()):
        S = [c2f_vertex[v] for v in cell if v in c2f_vertex]
        if len(S) >= 3:
            Sset = set(S)
            hits = []
            # triangles with all 3 vertices inside S (what extract_points(adjacent_cells=False) keeps)
            if len(Sset) >= 3:
                import itertools
                for trip in itertools.combinations(sorted(Sset), 3):
                    hits.extend(tri_by_set.get(frozenset(trip), []))
            if hits:
                cell_to_face[ci] = sorted(hits)
    return cell_to_face


def update_attributes(tet, cells, cell_to_face, graph, compute_in_air=True, dijkstra_cache=None):
    """reference.update_tet_attributes, minus the unused path_to_bottom array. Returns dijkstra results."""
    surface_mesh = tet.extract_surface(algorithm='dataset_surface')
    n = tet.number_of_cells
    fn = np.asarray(surface_mesh.face_normals)
    fc = np.asarray(surface_mesh.cell_centers().points)

    face_normal = np.full((n, 3), np.nan)
    face_center = np.full((n, 3), np.nan)
    for ci, fids in cell_to_face.items():
        nrm = fn[fids]
        face_normal[ci] = nrm[np.argmin(nrm[:, 2])]
        ctr = fc[fids]
        face_center[ci] = ctr[np.argmin(ctr[:, 2])]
    face_normal = face_normal / np.linalg.norm(face_normal, axis=1)[:, None]
    tet.cell_data['face_normal'] = face_normal
    tet.cell_data['face_center'] = face_center
    tet.cell_data["cell_center"] = tet.cell_centers().points
    cell_center = np.asarray(tet.cell_data["cell_center"])

    bottom_cell_threshold = np.nanmin(face_center[:, 2]) + 0.3
    with np.errstate(invalid='ignore'):
        bottom_cells_mask = face_center[:, 2] < bottom_cell_threshold
    tet.cell_data['is_bottom'] = bottom_cells_mask
    bottom_cells = np.where(bottom_cells_mask)[0]

    fnorm = face_normal.copy()
    fnorm[bottom_cells_mask] = np.nan
    overhang_angle = np.arccos(np.dot(fnorm, up_vector))
    tet.cell_data['overhang_angle'] = overhang_angle

    if dijkstra_cache is not None and dijkstra_cache[0] == set(bottom_cells.tolist()):
        distances, paths = dijkstra_cache[1], dijkstra_cache[2]  # same graph + same sources -> same result
    else:
        with TIMER("dijkstra (networkx)"):
            distances, paths = nx.multi_source_dijkstra(graph, set(bottom_cells))

    in_air = np.full(n, False)
    if compute_in_air:
        IN_AIR_THRESHOLD = 1
        z = cell_center[:, 2]
        # max z over each path: path[c] = path[parent] + [c] in networkx, so fold along parents
        maxz = np.full(n, -np.inf)
        order = sorted(paths, key=lambda c: len(paths[c]))
        for c in order:
            p = paths[c]
            if len(p) == 1:
                maxz[c] = z[c]
            else:
                parent = p[-2]
                zc = z[c]
                mp = maxz[parent]
                maxz[c] = zc if zc > mp else mp
        lens = np.array([len(paths[c]) for c in range(n)])
        in_air = (lens > 1) & (maxz > z + IN_AIR_THRESHOLD)
    tet.cell_data['in_air'] = in_air
    return bottom_cells_mask, bottom_cells, distances, paths


def prepare_mesh(model_path, part_offset=(0., 0., 0.)):
    """Notebook cell 2 (fast)."""
    ctx = MeshContext()
    with TIMER("load + tetgen"):
        tet = meshio_s4.load_and_tetrahedralize(model_path, part_offset)
    cells = meshio_s4.cells_of(tet)
    ctx.tet, ctx.cells = tet, cells

    with TIMER("neighbours (numpy + ordered sets)"):
        pairs = compute_neighbours(cells, tet.number_of_points)
        for kind, arr in pairs.items():
            tet.field_data[f"cell_{kind}_neighbours"] = arr
        ctx.pairs = pairs
        ctx.nbr = neighbour_dicts(pairs, tet.number_of_cells)

    with TIMER("networkx graph (bulk)"):
        cc = np.asarray(tet.cell_centers().points)
        pp = pairs["point"]
        diff = cc[pp[:, 0]] - cc[pp[:, 1]]
        dist = np.array([np.linalg.norm(d) for d in diff])  # same 1-D norm as the notebook
        g = nx.Graph()
        g.add_weighted_edges_from(zip(pp[:, 0].tolist(), pp[:, 1].tolist(), dist.tolist()))
        ctx.graph = g

    with TIMER("tet attributes (hashing)"):
        surface_mesh = tet.extract_surface(algorithm='dataset_surface')
        ctx.cell_to_face = cell_face_maps(tet, cells, surface_mesh)
        tet.add_field_data(base64.b64encode(pickle.dumps(ctx.cell_to_face)).decode('utf-8'), "cell_to_face")
        tet.add_field_data(cells, "cells")
        tet.add_field_data(tet.points, "cell_vertices")
        mask, bottom, dist_b, paths_b = update_attributes(tet, cells, ctx.cell_to_face, ctx.graph)
        tet.cell_data['overhang_angle'][bottom] = np.nan
        ctx.bottom_cells_mask, ctx.bottom_cells = mask, bottom
        ctx.distances_to_bottom, ctx.paths_to_bottom = dist_b, paths_b
    return ctx


# --------------------------------------------------------------------- rotation field

def planeFit(points):
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1))
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    M = np.dot(x, x.T)
    return ctr, svd(M)[0][:, -1]


def path_length_to_base_gradient(ctx, tet, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO):
    # tet: the mesh being deformed this iteration; ctx: original mesh (the notebook's globals:
    # neighbour graph, bottom_cells and its Dijkstra search are always those of the input mesh)
    n = tet.number_of_cells
    grad = np.zeros(n)
    cell_center = np.asarray(tet.cell_data["cell_center"])
    face_normal = np.asarray(tet.cell_data["face_normal"])
    nbr_edge = ctx.nbr["edge"]
    nbr_point = ctx.nbr["point"]

    with np.errstate(invalid='ignore'):
        is_overhang = np.arccos(face_normal @ np.array([0, 0, 1])) > np.deg2rad(90 + MAX_OVERHANG)
    candidates = np.nonzero(is_overhang & ~ctx.bottom_cells_mask)[0]
    dist = np.full(n, np.nan)
    closest = np.zeros(n, dtype=int)
    for c in candidates.tolist():
        closest[c] = ctx.paths_to_bottom[c][0]
        dist[c] = ctx.distances_to_bottom[c]
    tet.cell_data["cell_distance_to_bottom"] = dist

    with TIMER("gradient: plane fits"):
        for c in candidates.tolist():  # cells with non-nan distance
            local = np.hstack((nbr_edge[c], c))
            lpl = dist[local.astype(int)]
            valid = ~np.isnan(lpl)
            local = np.array(local)[valid]
            lpl = lpl[valid]
            if len(lpl) < 3:
                direction = cell_center[closest[c], :2] - cell_center[c, :2]
                direction /= np.linalg.norm(direction)
                ccn = cell_center[c, :2].copy()
                ccn /= np.linalg.norm(ccn)
                with np.errstate(invalid='ignore', divide='ignore'):
                    ord_ = np.dot(ccn, direction) / np.abs(np.dot(ccn, direction))
                if np.isnan(ord_):
                    ord_ = 0
                grad[c] = ord_
            else:
                local = local.astype(int)
                pts = np.hstack((cell_center[local, :2], lpl[:, None]))
                _, normal = planeFit(pts.T)
                ccn = cell_center[c, :2] / np.linalg.norm(cell_center[c, :2])
                gr = np.dot(ccn, normal[:2])
                if np.isnan(gr):
                    gr = np.mean(grad[local][~np.isnan(grad[local])])
                    if np.isnan(gr):
                        gr = 0
                grad[c] = gr

    with TIMER("gradient: smoothing (1 pass == notebook's 30)"):
        if INITIAL_ROTATION_FIELD_SMOOTHING != 0:
            smoothed = np.zeros(n)
            for c in np.nonzero(grad != 0)[0].tolist():
                neighbours = nbr_point[c]
                local = neighbours.copy()
                for nb in neighbours:
                    local.extend(nbr_point[nb])
                local = np.array(list(set(local)))
                local = local[grad[local] != 0]
                smoothed[c] = np.mean(grad[local])
            grad = smoothed

    if not SET_INITIAL_ROTATION_TO_ZERO:
        grad[grad == 0] = np.nan
    tet.cell_data["path_length_to_base_gradient"] = grad
    return grad


def initial_rotation_field(ctx, tet, p):
    irf = np.abs(np.deg2rad(90 + p["MAX_OVERHANG"]) - tet.cell_data['overhang_angle'])
    grad = path_length_to_base_gradient(ctx, tet, p["MAX_OVERHANG"], p["INITIAL_ROTATION_FIELD_SMOOTHING"], p["SET_INITIAL_ROTATION_TO_ZERO"])
    if p["STEEP_OVERHANG_COMPENSATION"]:
        ia = tet.cell_data["in_air"]
        irf[ia] += 2 * (np.deg2rad(180) - tet.cell_data['overhang_angle'][ia])
    irf *= grad
    irf = np.clip(irf * p["ROTATION_MULTIPLIER"], -np.deg2rad(360), np.deg2rad(360))
    irf = np.clip(irf, p["MAX_NEG_ROTATION"], p["MAX_POS_ROTATION"])
    tet.cell_data["initial_rotation_field"] = irf
    return irf


def optimize_rotations(ctx, tet, p, verbose=0):
    with TIMER("initial rotation field"):
        irf = initial_rotation_field(ctx, tet, p)
    valid = np.where(~np.isnan(irf))[0]
    W = p["NEIGHBOUR_LOSS_WEIGHT"]
    cfn = tet.field_data["cell_face_neighbours"]
    c1, c2 = cfn[:, 0], cfn[:, 1]
    E = len(cfn)
    n_rows = E + len(valid)
    n = tet.number_of_cells
    rows = np.concatenate([np.arange(E), np.arange(E), E + np.arange(len(valid))])
    cols = np.concatenate([c1, c2, valid])
    irf_valid = irf[valid]

    pattern = FastJacobian(n_rows, n, rows, cols)

    def fun(x):
        d = x[c1] - x[c2]
        return np.concatenate((W * d**2, (x[valid] - irf_valid)**2))

    def jac(x):
        d = x[c1] - x[c2]
        return pattern.build(np.concatenate([2 * W * d, -2 * W * d, 2 * (x[valid] - irf_valid)]))

    sparsity = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_rows, n))
    with TIMER("rotation least_squares (TRF)"), fast_trf(pattern):
        res = least_squares(fun, np.zeros(n), jac=jac, max_nfev=p["ROTATION_ITERATIONS"],
                            jac_sparsity=sparsity, verbose=verbose, method='trf', ftol=1e-6)
    return res.x


def rotation_matrices(tet, rotation_field):
    from scipy.spatial.transform import Rotation as R
    cxy = tet.cell_data["cell_center"][:, :2]
    c3 = np.hstack([cxy, np.zeros((cxy.shape[0], 1))])
    t = np.cross(np.array([0, 0, 1]), c3)
    with np.errstate(invalid='ignore'):
        t /= np.linalg.norm(t, axis=1)[:, None]
    t[np.isnan(t).any(axis=1)] = [1, 0, 0]
    return R.from_rotvec(rotation_field[:, None] * t).as_matrix()


def calculate_deformation(tet, rotation_field, iterations, verbose=0):
    cells = tet.field_data["cells"]
    n_cells = tet.number_of_cells
    n_pts = tet.number_of_points
    params = tet.points.copy().flatten()
    R_ = rotation_matrices(tet, rotation_field)
    old_vertices = tet.field_data["cell_vertices"][cells]
    old_T = np.einsum('ijk,ikl->ijl', R_, (N @ old_vertices).transpose(0, 2, 1))

    # Jacobian pattern: row = cell, cols = vertex*3+dim, in the order lil assignment produced
    cell_idx = np.repeat(np.arange(n_cells), cells.shape[1])
    vert_idx = np.ravel(cells)
    rows = np.concatenate([cell_idx] * 3)
    cols = np.concatenate([vert_idx * 3 + dim for dim in range(3)])

    pattern = FastJacobian(n_cells, len(params), rows, cols)
    cache = {"x": None, "diff": None}

    def _diff(x):
        # fun and jac are evaluated at the same x by TRF; compute (N V - R N V0) once
        if cache["x"] is not None and np.array_equal(cache["x"], x):
            return cache["diff"]
        nv = x[:n_pts * 3].reshape(-1, 3)
        nt = (N @ nv[cells]).transpose(0, 2, 1)
        cache["x"], cache["diff"] = x.copy(), nt - old_T
        return cache["diff"]

    def fun(x):
        return np.linalg.norm(_diff(x), axis=(1, 2))**2

    def jac(x):
        diff = _diff(x).transpose(0, 2, 1)
        return pattern.build(np.concatenate([(2 * diff[:, :, dim]).ravel() for dim in range(3)]))

    sparsity = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_cells, len(params)))
    with TIMER("deformation least_squares (TRF)"), fast_trf(pattern):
        res = least_squares(fun, params, max_nfev=iterations, verbose=verbose, jac=jac,
                            jac_sparsity=sparsity, method='trf', x_scale='jac')
    return res.x[:n_pts * 3].reshape(-1, 3)


def next_iteration_mesh(ctx, prev_tet, new_vertices):
    """Notebook cell 7 tail + cell 9: deformed mesh with attributes recomputed, ready for another iteration."""
    t = pv.UnstructuredGrid(ctx.tet.cells, np.full(ctx.tet.number_of_cells, pv.CellType.TETRA), new_vertices)
    t.field_data["cells"] = ctx.cells
    t.field_data["cell_vertices"] = t.points
    t.field_data["cell_face_neighbours"] = ctx.tet.field_data["cell_face_neighbours"]
    cached = (set(ctx.bottom_cells.tolist()), ctx.distances_to_bottom, ctx.paths_to_bottom)
    update_attributes(t, ctx.cells, ctx.cell_to_face, ctx.graph, dijkstra_cache=cached)
    return t


def deform(model_path, params=None, verbose=0):
    """Cells 2, (4, 7, [9]) per iteration, 11. Returns (ctx, deformed_grid_with_offset_applied, rotation_field)."""
    from .params import expand_iterations
    part_offset, iterations = expand_iterations(params)
    with TIMER("mesh setup"):
        ctx = prepare_mesh(model_path, part_offset)
    tet = ctx.tet
    for it, p in enumerate(iterations):
        tag = f" [iteration {it + 1}/{len(iterations)}]" if len(iterations) > 1 else ""
        with TIMER("optimize_rotations" + tag):
            rf = optimize_rotations(ctx, tet, p, verbose)
        with TIMER("calculate_deformation" + tag):
            nv = calculate_deformation(tet, rf, p["DEFORMATION_ITERATIONS"], verbose)
        if it < len(iterations) - 1:
            with TIMER("attributes for next iteration"):
                tet = next_iteration_mesh(ctx, tet, nv)
    d = pv.UnstructuredGrid(ctx.tet.cells, np.full(ctx.tet.number_of_cells, pv.CellType.TETRA), nv)
    x_min, x_max, y_min, y_max, z_min, z_max = d.bounds
    d.points -= np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
    return ctx, d, rf
