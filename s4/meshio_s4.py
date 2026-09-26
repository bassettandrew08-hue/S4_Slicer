"""Mesh loading and VTK helpers shared by the fast pipeline."""
import numpy as np
import open3d as o3d
import pyvista as pv
import tetgen


def load_and_tetrahedralize(model_path, part_offset=(0., 0., 0.)):
    """Notebook cell 2 (first part): identical calls, so identical tetgen output."""
    mesh = o3d.io.read_triangle_mesh(model_path)
    tg = tetgen.TetGen(np.asarray(mesh.vertices), np.asarray(mesh.triangles))
    tg.tetrahedralize()
    grid = tg.grid
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
