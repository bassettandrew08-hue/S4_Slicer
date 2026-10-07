"""Mesh loading and VTK helpers shared by the fast pipeline."""
import os
import time

import numpy as np
import open3d as o3d
import pyvista as pv
import tetgen
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


def clean_surface(vertices, faces, snap=1e-3, min_height=1e-6, min_volume=1e-3):
    """Remove CAD-export junk that TetGen reports as self-intersections: merge vertices closer than `snap` mm, drop
    collapsed and sliver triangles (height under `min_height` mm) and duplicate faces, then drop loose pieces
    enclosing less than `min_volume` mm^3. Vertices only move onto a neighbour within `snap` mm (chains of close
    vertices merge into one)."""
    v = np.asarray(vertices, dtype=float)
    pairs = cKDTree(v).query_pairs(snap, output_type="ndarray")
    n = len(v)
    _, cluster = connected_components(coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)),
                                      directed=False)
    _, keep = np.unique(cluster, return_index=True)  # each cluster becomes its first vertex
    v = v[keep]
    f = cluster[np.asarray(faces)]
    f = f[(f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 0] != f[:, 2])]
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    longest = np.max([np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1), np.linalg.norm(a - c, axis=1)], 0)
    f = f[area2 / longest > min_height]
    _, first = np.unique(np.sort(f, axis=1), axis=0, return_index=True)
    f = f[np.sort(first)]
    # loose flat pieces (no volume) left over from the slivers. A piece is faces joined by edges that exactly two
    # faces share: junk can touch the part at a vertex or hang off one of its edges
    m = len(f)
    edges = np.sort(np.stack([f, np.roll(f, -1, axis=1)], axis=2).reshape(-1, 2), axis=1)
    _, edge_id, edge_count = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
    edge_id = edge_id.reshape(-1)
    face_of = np.repeat(np.arange(m), 3)
    manifold = edge_count[edge_id] == 2
    adj = coo_matrix((np.ones(manifold.sum()), (face_of[manifold], m + edge_id[manifold])), shape=(4 * m, 4 * m))
    _, comp = connected_components(adj, directed=False)
    fc = comp[:m]
    vol = np.bincount(fc, np.einsum("ij,ij->i", v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]])) / 6)
    f = f[np.abs(vol[fc]) > min_volume]
    used, f = np.unique(f, return_inverse=True)
    return v[used], f.reshape(-1, 3)


def tetrahedralize(vertices, faces):
    """TetGen with default switches. If TetGen rejects the surface (self-intersections), retry once on the
    clean_surface version; meshes that tetrahedralize as they are never reach the cleanup."""
    t0 = time.time()
    try:
        tg = tetgen.TetGen(vertices, faces)
        tg.tetrahedralize()
    except RuntimeError as e:
        for fn in ("_skipped.face", "_skipped.node"):  # TetGen's dump of the bad faces, written into the cwd
            if os.path.exists(fn) and os.path.getmtime(fn) >= t0 - 1:
                os.remove(fn)
        v, f = clean_surface(vertices, faces)
        print(f"[mesh] TetGen: {e} Retrying after removing degenerate triangles "
              f"({len(faces)} -> {len(f)} faces)", flush=True)
        tg = tetgen.TetGen(v, f)
        tg.tetrahedralize()
    return tg.grid


def load_and_tetrahedralize(model_path, part_offset=(0., 0., 0.)):
    """Notebook cell 2 (first part): identical calls, so identical tetgen output."""
    mesh = o3d.io.read_triangle_mesh(model_path)
    grid = tetrahedralize(np.asarray(mesh.vertices), np.asarray(mesh.triangles))
    x_min, x_max, y_min, y_max, z_min, z_max = grid.bounds
    grid.points -= np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min]) + np.asarray(part_offset, dtype=float)
    return grid


def cells_of(grid):
    return grid.cells.reshape(-1, 5)[:, 1:]


def make_find_cells(deformed_grid):
    """Containing cell per point, falling back to the closest cell (as the notebook does)."""
    def find(P):
        containing = np.asarray(deformed_grid.find_containing_cell(P))
        miss = containing == -1
        if miss.any():
            containing = containing.copy()
            containing[miss] = np.asarray(deformed_grid.find_closest_cell(P[miss]))
        return containing
    return find


def deformed_grid_from_points(input_grid, new_points):
    return pv.UnstructuredGrid(input_grid.cells, np.full(input_grid.number_of_cells, pv.CellType.TETRA), new_points)
