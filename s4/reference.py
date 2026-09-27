"""
Faithful, un-optimized port of main.ipynb (Joshua Bird's S4 Slicer).

This module is the *golden reference* that every optimization is checked
against. The numerical code below is copied from the notebook cells with only
these mechanical changes:
  * notebook globals (cell_neighbour_graph, bottom_cells, cell_neighbour_dict,
    model_name, input_tet, ...) live in module globals set by the run_* helpers
  * plotting / .plot() calls removed; GIF writing kept behind SAVE_GIF
  * StageTimer blocks wrapped around the existing code (no logic changes)
  * the user's fixes are kept: MAX_POS/NEG_ROTATION = +/-360 deg, 3D zero-padded
    vectors for np.cross, no Linux-only calls

Do NOT optimize this file. Put faster code in s4/fast.py.
"""
import base64
import pickle
import time

import networkx as nx
import numpy as np
import open3d as o3d
import pyvista as pv
import tetgen
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix
from scipy.spatial.transform import Rotation as R

from .timing import TIMER

up_vector = np.array([0, 0, 1])

# notebook globals
model_name = None
cell_neighbour_graph = None
cell_neighbour_dict = None
bottom_cells = None
bottom_cells_mask = None
bottom_cell_groups = None
input_tet = None
undeformed_tet = None
save_gif_i = 0


def encode_object(obj):
    return base64.b64encode(pickle.dumps(obj)).decode('utf-8')


def decode_object(encoded_str):
    return pickle.loads(base64.b64decode(encoded_str))


# ----------------------------------------------------------------------------
# Cell 2
# ----------------------------------------------------------------------------

def update_tet_attributes(tet):
    '''
    Calculate face normals, face centers, cell centers, and overhang angles for each cell in the tetrahedral mesh.
    '''

    surface_mesh = tet.extract_surface(algorithm='dataset_surface')
    cell_to_face = decode_object(tet.field_data["cell_to_face"])

    # put general data in field_data for easy access
    cells = tet.cells.reshape(-1, 5)[:, 1:]  # assume all cells have 4 vertices
    tet.add_field_data(cells, "cells")
    cell_vertices = tet.points
    tet.add_field_data(cell_vertices, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]  # assume all faces have 3 vertices
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    tet.cell_data['face_normal'] = np.full((tet.number_of_cells, 3), np.nan)
    surface_mesh_face_normals = surface_mesh.face_normals
    for cell_index, face_indices in cell_to_face.items():
        face_normals = surface_mesh_face_normals[face_indices]
        # get the normal facing the most down
        most_down_normal_index = np.argmin(face_normals[:, 2])
        tet.cell_data['face_normal'][cell_index] = face_normals[most_down_normal_index]
    tet.cell_data['face_normal'] = tet.cell_data['face_normal'] / np.linalg.norm(tet.cell_data['face_normal'], axis=1)[:, None]

    tet.cell_data['face_center'] = np.empty((tet.number_of_cells, 3))
    tet.cell_data['face_center'][:, :] = np.nan
    surface_mesh_cell_centers = surface_mesh.cell_centers().points
    for cell_index, face_indices in cell_to_face.items():
        face_centers = surface_mesh_cell_centers[face_indices]
        # get the normal facing the most down
        most_down_center_index = np.argmin(face_centers[:, 2])
        tet.cell_data['face_center'][cell_index] = face_centers[most_down_center_index]

    tet.cell_data["cell_center"] = tet.cell_centers().points

    # calculate bottom cells
    bottom_cell_threshold = np.nanmin(tet.cell_data['face_center'][:, 2]) + 0.3
    bottom_cells_mask = tet.cell_data['face_center'][:, 2] < bottom_cell_threshold
    tet.cell_data['is_bottom'] = bottom_cells_mask
    bottom_cells = np.where(bottom_cells_mask)[0]

    face_normals = tet.cell_data['face_normal'].copy()
    face_normals[bottom_cells_mask] = np.nan  # make bottom faces not angled
    overhang_angle = np.arccos(np.dot(face_normals, up_vector))
    tet.cell_data['overhang_angle'] = overhang_angle

    overhang_direction = face_normals[:, :2].copy()
    with np.errstate(invalid='ignore'):
        overhang_direction /= np.linalg.norm(overhang_direction, axis=1)[:, None]
    tet.cell_data['overhang_direction'] = overhang_direction

    # calculate if cell will print in air by seeing if any cell centers along path to base are higher
    IN_AIR_THRESHOLD = 1
    tet.cell_data['in_air'] = np.full(tet.number_of_cells, False)

    with TIMER("update_tet_attributes: dijkstra"):
        _, paths_to_bottom = nx.multi_source_dijkstra(cell_neighbour_graph, set(bottom_cells))

    with TIMER("update_tet_attributes: path_to_bottom + in_air loops"):
        # put it in cell data
        tet.cell_data['path_to_bottom'] = np.full((tet.number_of_cells, np.max([len(x) for x in paths_to_bottom.values()])), -1)
        for cell_index, path_to_bottom in paths_to_bottom.items():
            tet.cell_data['path_to_bottom'][cell_index, :len(path_to_bottom)] = path_to_bottom

        # calculate if cell is in air
        for cell_index in range(tet.number_of_cells):
            path_to_bottom = paths_to_bottom[cell_index]
            if len(path_to_bottom) > 1:
                cell_heights = tet.cell_data['cell_center'][path_to_bottom, 2]
                if np.any(cell_heights > tet.cell_data['cell_center'][cell_index, 2] + IN_AIR_THRESHOLD):
                    tet.cell_data['in_air'][cell_index] = True

    return tet


def calculate_tet_attributes(tet):
    '''
    Calculate shared vertices between cells, cell to face & face to cell relations, and bottom cells of the tetrahedral mesh.
    '''

    surface_mesh = tet.extract_surface(algorithm='dataset_surface')

    # put general data in field_data for easy access
    cells = tet.cells.reshape(-1, 5)[:, 1:]  # assume all cells have 4 vertices
    tet.add_field_data(cells, "cells")
    cell_vertices = tet.points
    tet.add_field_data(cell_vertices, "cell_vertices")
    faces = surface_mesh.faces.reshape(-1, 4)[:, 1:]  # assume all faces have 3 vertices
    tet.add_field_data(faces, "faces")
    face_vertices = surface_mesh.points
    tet.add_field_data(face_vertices, "face_vertices")

    with TIMER("calculate_tet_attributes: shared_vertices loop (unused)"):
        # calculate shared vertices
        shared_vertices = []
        for cell_1, cell_2 in tet.field_data["cell_point_neighbours"]:
            shared_vertices_these_faces = np.intersect1d(cells[cell_1], cells[cell_2])
            for vertex in shared_vertices_these_faces:
                shared_vertices.append({
                    "cell_1_index": cell_1,
                    "cell_2_index": cell_2,
                    "cell_1_vertex_index": np.where(cells[cell_1] == vertex)[0][0],
                    "cell_2_vertex_index": np.where(cells[cell_2] == vertex)[0][0],
                })

    with TIMER("calculate_tet_attributes: cell<->face maps"):
        # calculate cell to face & face to cell relations
        cell_to_face = {}
        face_to_cell = {face_index: [] for face_index in range(len(faces))}
        cell_to_face_vertices = {}
        face_to_cell_vertices = {}
        for cell_vertex_index, cell_vertex in enumerate(tet.field_data["cell_vertices"].reshape(-1, 3)):
            face_vertex_index = np.where((face_vertices == cell_vertex).all(axis=1))[0]
            if len(face_vertex_index) == 1:
                cell_to_face_vertices[cell_vertex_index] = face_vertex_index[0]
                face_to_cell_vertices[face_vertex_index[0]] = cell_vertex_index

        for cell_index, cell in enumerate(tet.field_data["cells"]):
            face_vertex_indices = [cell_to_face_vertices[cell_vertex_index] for cell_vertex_index in cell if cell_vertex_index in cell_to_face_vertices]
            if len(face_vertex_indices) >= 3:
                extracted = surface_mesh.extract_points(face_vertex_indices, adjacent_cells=False)
                if extracted.number_of_cells >= 1:
                    cell_to_face[cell_index] = list(extracted.cell_data['vtkOriginalCellIds'])
                    for face_index in extracted.cell_data['vtkOriginalCellIds']:
                        face_to_cell[face_index].append(cell_index)

    tet.add_field_data(encode_object(cell_to_face), "cell_to_face")
    tet.add_field_data(encode_object(face_to_cell), "face_to_cell")

    # calculate has_face attribute
    tet.cell_data['has_face'] = np.zeros(tet.number_of_cells)
    for cell_index, face_indices in cell_to_face.items():
        tet.cell_data['has_face'][cell_index] = 1

    with TIMER("calculate_tet_attributes: update_tet_attributes"):
        tet = update_tet_attributes(tet)

    # calculate bottom cells
    bottom_cells_mask = tet.cell_data['is_bottom']
    bottom_cells = np.where(bottom_cells_mask)[0]

    tet.cell_data['overhang_angle'][bottom_cells] = np.nan

    return tet, bottom_cells_mask, bottom_cells


def run_cell2(model_path, PART_OFFSET=np.array([0., 0., 0.]), name=None):
    """Cell 2: load mesh, tetrahedralize, neighbours, tet attributes."""
    global model_name, cell_neighbour_graph, cell_neighbour_dict, bottom_cells, bottom_cells_mask
    global bottom_cell_groups, input_tet, undeformed_tet
    model_name = name

    with TIMER("load + tetgen"):
        mesh = o3d.io.read_triangle_mesh(model_path)
        # convert to tetrahedral mesh
        input_tet = tetgen.TetGen(np.asarray(mesh.vertices), np.asarray(mesh.triangles))
        input_tet.tetrahedralize()
        input_tet = input_tet.grid

    x_min, x_max, y_min, y_max, z_min, z_max = input_tet.bounds
    input_tet.points -= np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min]) + PART_OFFSET

    with TIMER("find neighbours (VTK cell_neighbors x3)"):
        # find neighbours
        cell_neighbour_dict = {neighbour_type: {face: [] for face in range(input_tet.number_of_cells)} for neighbour_type in ["point", "edge", "face"]}
        for neighbour_type in ["point", "edge", "face"]:
            cell_neighbours = []
            for cell_index in range(input_tet.number_of_cells):
                neighbours = input_tet.cell_neighbors(cell_index, f"{neighbour_type}s")
                for neighbour in neighbours:
                    if neighbour > cell_index:
                        cell_neighbours.append((cell_index, neighbour))
            for face_1, face_2 in np.array(cell_neighbours):
                cell_neighbour_dict[neighbour_type][face_1].append(face_2)
                cell_neighbour_dict[neighbour_type][face_2].append(face_1)

            input_tet.field_data[f"cell_{neighbour_type}_neighbours"] = np.array(cell_neighbours)

    with TIMER("build networkx neighbour graph"):
        cell_neighbour_graph = nx.Graph()
        cell_centers = input_tet.cell_centers().points
        for edge in input_tet.field_data["cell_point_neighbours"]:  # use point neighbours for best accuracy
            distance = np.linalg.norm(cell_centers[edge[0]] - cell_centers[edge[1]])
            cell_neighbour_graph.add_weighted_edges_from([(edge[0], edge[1], distance)])

    with TIMER("calculate_tet_attributes(input)"):
        bottom_cells_mask = None
        bottom_cells = None
        input_tet, bottom_cells_mask, bottom_cells = calculate_tet_attributes(input_tet)

    with TIMER("bottom cell groups"):
        # find bottom cell groups that are connected
        bottom_cell_graph = nx.Graph()
        for cell_index in bottom_cells:
            bottom_cell_graph.add_node(cell_index)
        cell_point_neighbour_dict = cell_neighbour_dict["point"]
        for cell_index in bottom_cells:
            for neighbour in cell_point_neighbour_dict[cell_index]:
                if neighbour in bottom_cells:
                    bottom_cell_graph.add_edge(cell_index, neighbour)

        bottom_cell_groups = [list(x) for x in list(nx.connected_components(bottom_cell_graph))]

    undeformed_tet = input_tet.copy()
    return input_tet


# ----------------------------------------------------------------------------
# Cell 3
# ----------------------------------------------------------------------------

def planeFit(points):
    """
    p, n = planeFit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    import numpy as np
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1))
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:, np.newaxis]
    M = np.dot(x, x.T)
    return ctr, svd(M)[0][:, -1]


def calculate_path_length_to_base_gradient(tet, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO):
    '''
    Calculate the path length to base gradient for each cell in the tetrahedral mesh with respect to the radial direction. This is used to determine the optimal rotation direction for each cell.

    returns: path_length_to_base_gradient. A scalar for each cell in the tetrahedral mesh. This is the gradient in the radial direction of the path length to the closest bottom cell.
    '''

    # calculate initial rotation direction for each face
    path_length_to_base_gradient = np.zeros((tet.number_of_cells))  # this is a scalar with respect to the radial direction. ie the vector pointing to the cell center

    with TIMER("gradient: dijkstra + overhang loop"):
        # find the path length for every overhang cell to a bottom cell
        cell_distance_to_bottom = np.empty((tet.number_of_cells))
        cell_distance_to_bottom[:] = np.nan
        distances_to_bottom, paths_to_bottom = nx.multi_source_dijkstra(cell_neighbour_graph, set(bottom_cells))
        closest_bottom_cell_indices = np.zeros((tet.number_of_cells), dtype=int)
        for cell_index in range(tet.number_of_cells):
            face_normal = tet.cell_data["face_normal"][cell_index]

            cell_is_overhang = np.arccos(np.dot(face_normal, [0, 0, 1])) > np.deg2rad(90 + MAX_OVERHANG)
            if cell_is_overhang and cell_index not in bottom_cells:
                closest_bottom_cell_indices[cell_index] = paths_to_bottom[cell_index][0]
                cell_distance_to_bottom[cell_index] = distances_to_bottom[cell_index]

        tet.cell_data["cell_distance_to_bottom"] = cell_distance_to_bottom

    with TIMER("gradient: per-cell plane fits"):
        # calculate the gradient of path length to base for each cell
        for cell_index in range(tet.number_of_cells):
            if not np.isnan(cell_distance_to_bottom[cell_index]):
                local_cells = cell_neighbour_dict["edge"][cell_index]
                local_cells = np.hstack((local_cells, cell_index))

                local_cell_path_lengths = [cell_distance_to_bottom[local_cell] for local_cell in local_cells]
                local_cell_path_lengths = np.array(local_cell_path_lengths)

                # remove neighbours with path length of nan
                local_cells = np.array(local_cells)[~np.isnan(local_cell_path_lengths)]
                local_cell_path_lengths = local_cell_path_lengths[~np.isnan(local_cell_path_lengths)]

                # if there are less than 3 neighbours with path length, roll to the closest bottom cell
                if len(local_cell_path_lengths) < 3:
                    location_to_roll_to = tet.cell_data["cell_center"][closest_bottom_cell_indices[cell_index], :2]

                    direction_to_bottom = location_to_roll_to - tet.cell_data["cell_center"][cell_index, :2]
                    direction_to_bottom /= np.linalg.norm(direction_to_bottom)

                    cell_center = tet.cell_data["cell_center"][cell_index, :2].copy()
                    cell_center /= np.linalg.norm(cell_center)

                    with np.errstate(invalid='ignore', divide='ignore'):
                        optimal_rotation_direction = np.dot(cell_center, direction_to_bottom) / np.abs(np.dot(cell_center, direction_to_bottom))
                    if np.isnan(optimal_rotation_direction):
                        optimal_rotation_direction = 0

                    path_length_to_base_gradient[cell_index] = optimal_rotation_direction

                # if there are 3 or more neighbours with path length, calculate the gradient in the radial direction
                # and use that as the optimal rotation direction
                else:
                    points = np.hstack((tet.cell_data["cell_center"][local_cells, :2], local_cell_path_lengths[:, None]))
                    _, plane_normal = planeFit(points.T)

                    cell_center_direction_normalized = tet.cell_data["cell_center"][cell_index, :2] / np.linalg.norm(tet.cell_data["cell_center"][cell_index, :2])
                    gradient_in_radial_direction = np.dot(cell_center_direction_normalized, plane_normal[:2])

                    # if the gradient is nan, use the average of the neighbours
                    if np.isnan(gradient_in_radial_direction):
                        gradient_in_radial_direction = np.mean(path_length_to_base_gradient[local_cells][~np.isnan(path_length_to_base_gradient[local_cells])])
                        if np.isnan(gradient_in_radial_direction):
                            gradient_in_radial_direction = 0

                    path_length_to_base_gradient[cell_index] = gradient_in_radial_direction

    with TIMER("gradient: smoothing loop"):
        # smooth path_length_to_base_gradient with neighbours
        # not needed because we do neighbour difference minimization in the optimization step?
        if INITIAL_ROTATION_FIELD_SMOOTHING != 0:
            for i in range(INITIAL_ROTATION_FIELD_SMOOTHING):
                smoothed_path_length_to_base_gradient = np.zeros((tet.number_of_cells))
                for cell_index in range(tet.number_of_cells):
                    if path_length_to_base_gradient[cell_index] != 0:
                        neighbours = cell_neighbour_dict["point"][cell_index]
                        local_cells = neighbours.copy()
                        for neighbour in neighbours:
                            local_cells.extend(cell_neighbour_dict["point"][neighbour])
                        local_cells = np.array(list(set(local_cells)))
                        local_cells = local_cells[path_length_to_base_gradient[local_cells] != 0]
                        smoothed_path_length_to_base_gradient[cell_index] = np.mean(path_length_to_base_gradient[local_cells])

            path_length_to_base_gradient = smoothed_path_length_to_base_gradient

    # replace 0 with nan
    if not SET_INITIAL_ROTATION_TO_ZERO:
        path_length_to_base_gradient[path_length_to_base_gradient == 0] = np.nan
    tet.cell_data["path_length_to_base_gradient"] = path_length_to_base_gradient  # very sexy

    return path_length_to_base_gradient


# ----------------------------------------------------------------------------
# Cell 4
# ----------------------------------------------------------------------------

def calculate_initial_rotation_field(tet, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION):
    '''
    Calculate the initial rotation field for each cell in the tetrahedral mesh to make overhangs less than MAX_OVERHANG.
    The direction of rotation ensures the part is printable.
    '''

    # create initial rotation field rotating faces to be in safe printing angle
    initial_rotation_field = np.full((tet.number_of_cells), np.nan)
    initial_rotation_field = np.abs(np.deg2rad(90 + MAX_OVERHANG) - tet.cell_data['overhang_angle'])

    path_length_to_base_gradient = calculate_path_length_to_base_gradient(tet, MAX_OVERHANG, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO)

    if STEEP_OVERHANG_COMPENSATION:
        initial_rotation_field[tet.cell_data["in_air"]] += 2 * (np.deg2rad(180) - tet.cell_data['overhang_angle'][tet.cell_data["in_air"]])

    # # Apply the path_length_to_base_gradient (optimal overhang direction) to the initial rotation field
    initial_rotation_field *= path_length_to_base_gradient

    # apply rotation multiplier
    initial_rotation_field = np.clip(initial_rotation_field * ROTATION_MULTIPLIER, -np.deg2rad(360), np.deg2rad(360))

    # clip to max rotation
    initial_rotation_field = np.clip(initial_rotation_field, MAX_NEG_ROTATION, MAX_POS_ROTATION)

    tet.cell_data["initial_rotation_field"] = initial_rotation_field

    return initial_rotation_field


def calculate_rotation_matrices(tet, rotation_field):
    '''
    Calculate the rotation matrices for each cell in the tetrahedral mesh given the scalar
    rotation field that gives a rotation for each cell. Cells are rotated around the axis
    perpendicular to the radial direction and the z-axis.
    '''

    # create rotation matrix from theta around axis
    cell_centers_xy = tet.cell_data["cell_center"][:, :2]
    cell_centers_3d = np.hstack([cell_centers_xy, np.zeros((cell_centers_xy.shape[0], 1))])
    tangential_vectors = np.cross(np.array([0, 0, 1]), cell_centers_3d)
    # normalize
    with np.errstate(invalid='ignore'):
        tangential_vectors /= np.linalg.norm(tangential_vectors, axis=1)[:, None]
    # replace nan with [1,0,0]
    tangential_vectors[np.isnan(tangential_vectors).any(axis=1)] = [1, 0, 0]

    rotation_matrices = R.from_rotvec(rotation_field[:, None] * tangential_vectors).as_matrix()

    return rotation_matrices


def calculate_unique_vertices_rotated(tet, rotation_field):
    rotation_matrices = calculate_rotation_matrices(tet, rotation_field)

    # rotate each face by the rotation field around its center
    unique_vertices = np.zeros((tet.number_of_cells, 4, 3))
    for cell_index, cell in enumerate(tet.field_data["cells"]):
        unique_vertices[cell_index] = tet.field_data["cell_vertices"][cell]

    cell_centers = tet.cell_data["cell_center"]

    unique_vertices_rotated = cell_centers.reshape(-1, 1, 3, 1) + rotation_matrices.reshape(-1, 1, 3, 3) @ (unique_vertices.reshape(-1, 4, 3, 1) - cell_centers.reshape(-1, 1, 3, 1))

    return unique_vertices_rotated


def apply_rotation_field_unique_vertices(tet, rotation_field):
    unique_vertices_rotated = calculate_unique_vertices_rotated(tet, rotation_field)

    unique_cells = np.zeros((tet.number_of_cells, 5), dtype=int)
    unique_cells[:, 0] = 4
    unique_cells[:, 1:] = np.arange(tet.number_of_cells * 4).reshape(-1, 4)

    new_tet = pv.UnstructuredGrid(unique_cells.flatten(), np.full(tet.number_of_cells, pv.CellType.TETRA), unique_vertices_rotated.reshape(-1, 3))

    return new_tet


def optimize_rotations(tet, NEIGHBOUR_LOSS_WEIGHT, MAX_OVERHANG, ROTATION_MULTIPLIER, ITERATIONS, SAVE_GIF, STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION, verbose=2):
    '''
    Optimize the rotation field for each cell in the tetrahedral mesh to make overhangs less
    than MAX_OVERHANG while keeping the rotation field smooth.
    '''

    plotter = None
    if SAVE_GIF:
        plotter = pv.Plotter(off_screen=True)
        plotter.open_gif(f'gifs/{model_name}_optimize_rotations.gif')

    with TIMER("initial rotation field"):
        initial_rotation_field = calculate_initial_rotation_field(tet, MAX_OVERHANG, ROTATION_MULTIPLIER, STEEP_OVERHANG_COMPENSATION, INITIAL_ROTATION_FIELD_SMOOTHING, SET_INITIAL_ROTATION_TO_ZERO, MAX_POS_ROTATION, MAX_NEG_ROTATION)
    num_cells_with_initial_rotation = np.sum(~np.isnan(initial_rotation_field))

    def save_gif(rotation_field):
        new_tet = apply_rotation_field_unique_vertices(tet, rotation_field)
        new_tet.cell_data["rotation_field"] = rotation_field
        mesh_actor = plotter.add_mesh(new_tet, clim=[-np.pi / 4, np.pi / 4], scalars="rotation_field", lighting=False)
        plotter.write_frame()
        plotter.remove_actor(mesh_actor)

    def objective_function(rotation_field):
        if SAVE_GIF:
            save_gif(rotation_field)

        cell_face_neighbours = tet.field_data["cell_face_neighbours"]
        neighbour_differences = rotation_field[cell_face_neighbours[:, 0]] - rotation_field[cell_face_neighbours[:, 1]]
        neighbour_losses = NEIGHBOUR_LOSS_WEIGHT * neighbour_differences**2

        valid_cell_indices = np.where(~np.isnan(initial_rotation_field))[0]
        initial_rotation_losses = (rotation_field[valid_cell_indices] - initial_rotation_field[valid_cell_indices])**2

        return np.concatenate((neighbour_losses, initial_rotation_losses))

    def objective_jacobian(rotation_field):
        cell_face_neighbours = tet.field_data["cell_face_neighbours"]
        jac = lil_matrix((len(cell_face_neighbours) + num_cells_with_initial_rotation, tet.number_of_cells), dtype=np.float32)

        cell_1 = cell_face_neighbours[:, 0]
        cell_2 = cell_face_neighbours[:, 1]

        differences = rotation_field[cell_1] - rotation_field[cell_2]

        jac[range(len(cell_face_neighbours)), cell_1] = 2 * NEIGHBOUR_LOSS_WEIGHT * differences
        jac[range(len(cell_face_neighbours)), cell_2] = -2 * NEIGHBOUR_LOSS_WEIGHT * differences

        valid_cell_indices = np.where(~np.isnan(initial_rotation_field))[0]

        jac[len(cell_face_neighbours) + np.arange(len(valid_cell_indices)), valid_cell_indices] = \
            2 * (rotation_field[valid_cell_indices] - initial_rotation_field[valid_cell_indices])

        return jac.tocsr()

    def jac_sparsity():
        cell_face_neighbours = tet.field_data["cell_face_neighbours"]
        sparsity = lil_matrix((len(cell_face_neighbours) + num_cells_with_initial_rotation, tet.number_of_cells), dtype=np.int8)

        for i, (cell_1, cell_2) in enumerate(cell_face_neighbours):
            sparsity[i, cell_1] = 1
            sparsity[i, cell_2] = 1

        valid_cell_indices = np.where(~np.isnan(initial_rotation_field))[0]
        i = 0
        for cell_index, initial_rotation in enumerate(initial_rotation_field):
            if cell_index in valid_cell_indices:
                sparsity[len(cell_face_neighbours) + i, cell_index] = 1
                i += 1

        return sparsity.tocsr()

    smoothed_rotation_field = np.zeros((tet.number_of_cells))

    with TIMER("rotation jac_sparsity build"):
        sparsity = jac_sparsity()
    with TIMER("rotation least_squares (TRF)"):
        result = least_squares(objective_function,
                               smoothed_rotation_field,
                               jac=objective_jacobian,
                               max_nfev=ITERATIONS,
                               jac_sparsity=sparsity,
                               verbose=verbose,
                               method='trf',
                               ftol=1e-6,
                               )

    if SAVE_GIF:
        plotter.close()

    return result.x


# ----------------------------------------------------------------------------
# Cell 7
# ----------------------------------------------------------------------------

N = np.eye(4) - 1 / 4 * np.ones((4, 4))  # the N matrix centers the vertices of a tetrahedron around the origin


def calculate_deformation(tet, rotation_field, ITERATIONS, SAVE_GIF, verbose=2):
    new_vertices = tet.points.copy()

    params = new_vertices.flatten()

    rotation_matrices = calculate_rotation_matrices(tet, rotation_field)

    old_vertices = tet.field_data["cell_vertices"][tet.field_data["cells"]]
    old_vertices_transformed = np.einsum('ijk,ikl->ijl', rotation_matrices, (N @ old_vertices).transpose(0, 2, 1))

    plotter = None
    if SAVE_GIF:
        plotter = pv.Plotter(off_screen=True)
        plotter.open_gif(f'gifs/{model_name}_calculate_deformation.gif')

    def save_gif(new_vertices):
        global save_gif_i
        save_gif_i += 1

        if save_gif_i % 10 != 0:
            return

        new_tet = pv.UnstructuredGrid(tet.cells, np.full(tet.number_of_cells, pv.CellType.TETRA), new_vertices)
        mesh_actor = plotter.add_mesh(new_tet)
        plotter.write_frame()
        plotter.remove_actor(mesh_actor)

    def objective_function(params):
        new_vertices = params[:tet.number_of_points * 3].reshape(-1, 3)

        if SAVE_GIF:
            save_gif(new_vertices)

        new_vertices_transformed = (N @ new_vertices[tet.field_data["cells"]]).transpose(0, 2, 1)

        position_losses = np.linalg.norm(new_vertices_transformed - old_vertices_transformed, axis=(1, 2))**2

        return position_losses

    def objective_jacobian(params):
        J = lil_matrix((tet.number_of_cells, len(params)), dtype=np.float32)

        new_vertices = params[:tet.number_of_points * 3].reshape(-1, 3)

        old_vertices = tet.field_data["cell_vertices"][tet.field_data["cells"]]

        new_vertices_transformed = (N @ new_vertices[tet.field_data["cells"]]).transpose(0, 2, 1)

        diff = new_vertices_transformed - old_vertices_transformed

        diff = diff.transpose(0, 2, 1)

        cell_indices = np.repeat(np.arange(tet.number_of_cells), len(tet.field_data["cells"][0]))
        vertex_indices = np.ravel(tet.field_data["cells"])

        for dim in range(3):
            J[cell_indices, vertex_indices * 3 + dim] = 2 * diff[:, :, dim].ravel()

        return J.tocsr()

    def jac_sparsity():
        sparsity = lil_matrix((tet.number_of_cells, len(params)), dtype=np.int8)

        cell_indices = np.repeat(np.arange(tet.number_of_cells), len(tet.field_data["cells"][0]))
        vertex_indices = np.ravel(tet.field_data["cells"])

        for dim in range(3):
            sparsity[cell_indices, vertex_indices * 3 + dim] = 1

        return sparsity.tocsr()

    with TIMER("deformation jac_sparsity build"):
        sparsity = jac_sparsity()
    with TIMER("deformation least_squares (TRF)"):
        result = least_squares(objective_function,
                               params,
                               max_nfev=ITERATIONS,
                               verbose=verbose,
                               jac=objective_jacobian,
                               jac_sparsity=sparsity,
                               method='trf',
                               x_scale='jac',
                               )

    if SAVE_GIF:
        plotter.close()

    return result.x[:tet.number_of_points * 3].reshape(-1, 3)


def deform(model_path, name, params, save_gif=False, verbose=2):
    """Cells 2, then (4, 7, [9]) per iteration, then 11. Returns (input_tet, deformed_tet, rotation_field)."""
    global undeformed_tet
    from .params import expand_iterations
    part_offset, iterations = expand_iterations(params)
    with TIMER("cell 2: mesh setup"):
        run_cell2(model_path, np.asarray(part_offset, dtype=float), name)

    for it, p in enumerate(iterations):
        tag = f" [iteration {it + 1}/{len(iterations)}]" if len(iterations) > 1 else ""
        with TIMER("cell 4: optimize_rotations" + tag):
            rotation_field = optimize_rotations(
                undeformed_tet,
                p["NEIGHBOUR_LOSS_WEIGHT"], p["MAX_OVERHANG"], p["ROTATION_MULTIPLIER"], p["ROTATION_ITERATIONS"],
                save_gif, p["STEEP_OVERHANG_COMPENSATION"], p["INITIAL_ROTATION_FIELD_SMOOTHING"],
                p["SET_INITIAL_ROTATION_TO_ZERO"], p["MAX_POS_ROTATION"], p["MAX_NEG_ROTATION"], verbose=verbose)

        with TIMER("cell 7: calculate_deformation" + tag):
            if p.get("DEFORMATION_METHOD", "island_free") == "notebook":
                new_vertices = calculate_deformation(undeformed_tet, rotation_field, p["DEFORMATION_ITERATIONS"], save_gif, verbose=verbose)
            else:  # not in the notebook: shared with the fast implementation (s4/island_free.py)
                from .fast_deform import deformation_step
                new_vertices = deformation_step(undeformed_tet, rotation_field, p, last=(it == len(iterations) - 1))
            deformed_tet = pv.UnstructuredGrid(undeformed_tet.cells, np.full(undeformed_tet.number_of_cells, pv.CellType.TETRA), new_vertices)

            for key in undeformed_tet.field_data.keys():
                deformed_tet.field_data[key] = undeformed_tet.field_data[key]
            for key in undeformed_tet.cell_data.keys():
                deformed_tet.cell_data[key] = undeformed_tet.cell_data[key]
            with TIMER("update_tet_attributes(deformed)"):
                deformed_tet = update_tet_attributes(deformed_tet)

        if it < len(iterations) - 1:
            undeformed_tet = deformed_tet.copy()  # cell 9

    # cell 11: make origin center bottom of bounding box
    x_min, x_max, y_min, y_max, z_min, z_max = deformed_tet.bounds
    offsets_applied = np.array([(x_min + x_max) / 2, (y_min + y_max) / 2, z_min])
    deformed_tet.points -= offsets_applied

    return input_tet, deformed_tet, rotation_field


def save_deformed_stl(deformed_tet, path):
    """Cell 11 (save part)."""
    deformed_tet.extract_surface(algorithm='dataset_surface').save(path)


# ----------------------------------------------------------------------------
# Cells 15-18
# ----------------------------------------------------------------------------

def tetrahedron_volume(p1, p2, p3, p4):
    mat = np.vstack([p2 - p1, p3 - p1, p4 - p1])
    return np.abs(np.linalg.det(mat)) / 6


def calc_barycentric_coordinates(tet_a, tet_b, tet_c, tet_d, point):
    total_volume = tetrahedron_volume(tet_a, tet_b, tet_c, tet_d)

    if total_volume == 0:
        raise ValueError("The points do not form a valid tetrahedron (zero volume).")

    vol_a = tetrahedron_volume(point, tet_b, tet_c, tet_d)
    vol_b = tetrahedron_volume(point, tet_a, tet_c, tet_d)
    vol_c = tetrahedron_volume(point, tet_a, tet_b, tet_d)
    vol_d = tetrahedron_volume(point, tet_a, tet_b, tet_c)

    lambda_a = vol_a / total_volume
    lambda_b = vol_b / total_volume
    lambda_c = vol_c / total_volume
    lambda_d = vol_d / total_volume

    return np.array([lambda_a, lambda_b, lambda_c, lambda_d])


def project_point_onto_plane(plane_x_axis, plane_y_axis, point):
    projected_x = np.sum(plane_x_axis * point, axis=1)
    projected_y = np.sum(plane_y_axis * point, axis=1)

    return np.array([projected_x, projected_y]).T


MAPPING_DEFAULTS = dict(
    SEG_SIZE=0.6,  # mm
    MAX_ROTATION=30,  # degrees
    MIN_ROTATION=-130,  # degrees
    NOZZLE_OFFSET=42,  # mm actuallt 41.5
    ROTATION_AVERAGING_ALPHA=0.2,
    RETRACTION_LENGTH=1.0,
    ROTATION_MAX_DELTA=np.deg2rad(1),
    MAX_EXTRUSION_MULTIPLIER=10,
    SPLIT_RETRACTIONS=False,  # not in the notebook; see map_gcode
    SMOOTH_EXTRUSION_MULTIPLIER=False,  # not in the notebook; see s4/fast_map.py
    EXTRUSION_MULTIPLIER_RANGE=None,  # not in the notebook; see s4/fast_map.py
    SAFE_TRAVEL_TRANSITIONS=False,  # not in the notebook; see s4/fast_map.py
)


def read_gcode_points(gcode_path, SEG_SIZE):
    """Cell 17: G-code reading/segmentation part (pygcode)."""
    from pygcode import Line
    pos = np.array([0., 0., 20.])
    feed = 5000
    gcode_points = []
    with open(gcode_path, 'r') as fh:
        for line_text in fh.readlines():
            line = Line(line_text)

            if not line.block.gcodes:
                continue

            for gcode in sorted(line.block.gcodes):
                if gcode.word == "G01" or gcode.word == "G00":
                    prev_pos = pos.copy()

                    if gcode.X is not None:
                        pos[0] = gcode.X
                    if gcode.Y is not None:
                        pos[1] = gcode.Y
                    if gcode.Z is not None:
                        pos[2] = gcode.Z

                    inv_time_feed = None
                    for word in line.block.words:
                        if word.letter == "F":
                            feed = word.value

                    extrusion = None
                    for param in line.block.modal_params:
                        if param.letter == "E":
                            extrusion = param.value

                    delta_pos = pos - prev_pos
                    distance = np.linalg.norm(delta_pos)
                    if distance > 0:
                        num_segments = -(-distance // SEG_SIZE)  # hacky round up
                        seg_distance = distance / num_segments

                        time_to_complete_move = (1 / feed) * seg_distance  # min/mm * mm = min
                        if time_to_complete_move == 0:
                            inv_time_feed = None
                        else:
                            inv_time_feed = 1 / time_to_complete_move  # 1/min

                        for i in range(int(num_segments)):
                            gcode_points.append({
                                "position": (prev_pos + delta_pos * (i + 1) / num_segments),
                                "command": gcode.word,
                                "extrusion": extrusion / num_segments if extrusion is not None else None,
                                "inv_time_feed": inv_time_feed,
                                "move_length": seg_distance,
                                "start_position": prev_pos,
                                "end_position": pos,
                                "unsegmented_move_length": distance,
                                "after_retract": False,
                                "feed": feed
                            })
                    else:
                        time_to_complete_move = (1 / feed) * distance  # min/mm * mm = min
                        if time_to_complete_move == 0:
                            inv_time_feed = None
                        else:
                            inv_time_feed = 1 / time_to_complete_move  # 1/min

                        gcode_points.append({
                            "position": pos.copy(),
                            "command": gcode.word,
                            "extrusion": extrusion,
                            "inv_time_feed": inv_time_feed,
                            "move_length": distance,
                            "unsegmented_move_length": distance,
                            "after_retract": False,
                            "feed": feed
                        })
    # the notebook's gcode.word is a pygcode Word; normalise to str for comparisons downstream
    for p in gcode_points:
        p["command"] = str(p["command"])
    return gcode_points


def map_gcode(input_tet, deformed_tet, gcode_path, mp=None):
    """Cells 16 + 17. Returns (new_gcode_points, stats)."""
    mp = {**MAPPING_DEFAULTS, **(mp or {})}
    SEG_SIZE = mp["SEG_SIZE"]
    MAX_ROTATION = mp["MAX_ROTATION"]
    MIN_ROTATION = mp["MIN_ROTATION"]

    with TIMER("cell 16: calculate_tet_attributes(deformed)"):
        deformed_tet, _, _ = calculate_tet_attributes(deformed_tet)

    with TIMER("per-cell rotations (2D Kabsch loop)"):
        vertex_transformations = deformed_tet.points - input_tet.points

        cell_centers_xy = input_tet.cell_data["cell_center"][:, :2]
        cell_centers_3d = np.hstack([cell_centers_xy, np.zeros((cell_centers_xy.shape[0], 1))])
        tangential_vectors = np.cross(np.array([0, 0, 1]), cell_centers_3d)
        with np.errstate(invalid='ignore'):
            tangential_vectors /= np.linalg.norm(tangential_vectors, axis=1)[:, None]
        tangential_vectors[np.isnan(tangential_vectors).any(axis=1)] = [1, 0, 0]

        num_cells_per_vertex = np.zeros((input_tet.number_of_points))
        for cell_index, cell in enumerate(input_tet.field_data["cells"]):
            num_cells_per_vertex[cell] += 1
        vertex_rotations = np.zeros((deformed_tet.number_of_points))
        cell_rotations = np.zeros((deformed_tet.number_of_cells))
        for cell_index, cell in enumerate(deformed_tet.field_data["cells"]):
            new_vertices = deformed_tet.field_data["cell_vertices"][cell]
            new_cell_center = deformed_tet.cell_data["cell_center"][cell_index]
            old_vertices = input_tet.field_data["cell_vertices"][cell]
            old_cell_center = input_tet.cell_data["cell_center"][cell_index]

            new_vertices -= new_cell_center
            old_vertices -= old_cell_center

            plane_x_vector = old_cell_center[:2] / np.linalg.norm(old_cell_center[:2])
            plane_x_vector = np.array([plane_x_vector[0], plane_x_vector[1], 0])
            plane_y_vector = np.array([0, 0, 1])

            new_vertices_projected = project_point_onto_plane(plane_x_vector, plane_y_vector, new_vertices)
            old_vertices_projected = project_point_onto_plane(plane_x_vector, plane_y_vector, old_vertices)

            covariance_matrix = np.dot(new_vertices_projected.T, old_vertices_projected)
            U, _, Vt = np.linalg.svd(covariance_matrix)
            rotation_matrix = np.dot(U, Vt)

            rotation = -np.arccos(min(max(rotation_matrix[0, 0], -1), 1))
            if rotation_matrix[1, 0] < 0:
                rotation = -rotation

            rotation = max(min(rotation, np.deg2rad(MAX_ROTATION)), np.deg2rad(MIN_ROTATION))

            cell_rotations[cell_index] = rotation

            for vertex_index in cell:
                vertex_rotations[vertex_index] += rotation / num_cells_per_vertex[vertex_index]

    with TIMER("z squish scales loop"):
        tet_rotation_matrices = calculate_rotation_matrices(input_tet, cell_rotations)
        z_squish_scales = np.full((deformed_tet.number_of_cells), np.nan)
        cell_vol0 = np.full((deformed_tet.number_of_cells), np.nan)  # not in the notebook (smoothed multiplier)
        cell_vold = np.full((deformed_tet.number_of_cells), np.nan)
        for cell_index, cell in enumerate(deformed_tet.field_data["cells"]):
            warped_vertices = deformed_tet.field_data["cell_vertices"][cell]
            unwarped_vertices = input_tet.field_data["cell_vertices"][cell]

            unwarped_vertices_rotated = (tet_rotation_matrices[cell_index].reshape(1, 3, 3) @ unwarped_vertices.reshape(4, 3, 1)).reshape(4, 3)

            z_squish_scales[cell_index] = tetrahedron_volume(*unwarped_vertices) / tetrahedron_volume(*warped_vertices)
            cell_vol0[cell_index] = tetrahedron_volume(*unwarped_vertices)
            cell_vold[cell_index] = tetrahedron_volume(*warped_vertices)
        # per-vertex volume ratio (same formula as fast_map.vertex_volume_ratio)
        _v0 = np.zeros(deformed_tet.number_of_points)
        _vd = np.zeros(deformed_tet.number_of_points)
        np.add.at(_v0, deformed_tet.field_data["cells"].ravel(), np.repeat(cell_vol0, 4))
        np.add.at(_vd, deformed_tet.field_data["cells"].ravel(), np.repeat(cell_vold, 4))
        vertex_volume_ratio = _v0 / _vd

    with TIMER("read + segment gcode (pygcode)"):
        gcode_points = read_gcode_points(gcode_path, SEG_SIZE)

    with TIMER("find_containing_cell + find_closest_cell"):
        gcode_points_containing_cells = deformed_tet.find_containing_cell([point["position"] for point in gcode_points])
        gcode_points_closest_cells = deformed_tet.find_closest_cell([point["position"] for point in gcode_points])
        gcode_points_containing_cells[gcode_points_containing_cells == -1] = gcode_points_closest_cells[gcode_points_containing_cells == -1]

    with TIMER("per-point transform loop (barycentric)"):
        new_gcode_points = []
        prev_new_position = None
        travelling_over_air = False
        travelling = False
        prev_position = None
        prev_rotation = 0
        prev_travelling = False
        prev_command = "G00"
        ROTATION_AVERAGING_ALPHA = mp["ROTATION_AVERAGING_ALPHA"]
        RETRACTION_LENGTH = mp["RETRACTION_LENGTH"]
        ROTATION_MAX_DELTA = mp["ROTATION_MAX_DELTA"]
        MAX_EXTRUSION_MULTIPLIER = mp["MAX_EXTRUSION_MULTIPLIER"]
        SPLIT_RETRACTIONS = mp["SPLIT_RETRACTIONS"]
        SMOOTH_EXTRUSION_MULTIPLIER = mp["SMOOTH_EXTRUSION_MULTIPLIER"]
        EXTRUSION_MULTIPLIER_RANGE = mp["EXTRUSION_MULTIPLIER_RANGE"]
        SAFE_TRAVEL_TRANSITIONS = mp["SAFE_TRAVEL_TRANSITIONS"]
        last_bary = {}
        lost_vertices = []
        highest_printed_point = 0
        no_cell_positions = []
        bad_barycentric_positions = []
        for cell_index, (gcode_point, containing_cell_index) in enumerate(zip(gcode_points, gcode_points_containing_cells)):
            position = gcode_point["position"]
            command = gcode_point["command"]
            inv_time_feed = gcode_point["inv_time_feed"]
            extrusion = gcode_point["extrusion"]

            def barycentric_interpolate_to_get_new_position_and_rotation(position, containing_cell_index, command, cell_index):
                if command == "G00" and containing_cell_index == -1:
                    no_cell_positions.append(position)
                    return None, None
                if command == "G01" and containing_cell_index == -1:
                    containing_cell_index = gcode_points_closest_cells[cell_index]

                vertiex_indices = deformed_tet.field_data["cells"][containing_cell_index]
                cell_vertices = deformed_tet.field_data["cell_vertices"][vertiex_indices]
                barycentric_coordinates = calc_barycentric_coordinates(cell_vertices[0], cell_vertices[1], cell_vertices[2], cell_vertices[3], position)

                if np.sum(barycentric_coordinates) > 1.01:
                    bad_barycentric_positions.append(position)
                    return None, None

                transformation = vertex_transformations[vertiex_indices] * barycentric_coordinates[:, None]
                transformation = np.sum(transformation, axis=0)
                new_position = position - transformation

                rotation = np.sum(vertex_rotations[vertiex_indices] * barycentric_coordinates)

                last_bary["v"] = (vertiex_indices, barycentric_coordinates)
                return new_position, rotation

            dont_smooth_rotation = False
            reentry = None
            new_position, rotation = barycentric_interpolate_to_get_new_position_and_rotation(position, containing_cell_index, command, cell_index)
            if new_position is None:
                if command == "G01":
                    lost_vertices.append(position)
                    continue
                elif command == "G00" and not travelling_over_air and prev_new_position is not None:
                    new_position = np.array([prev_new_position[0], prev_new_position[1], highest_printed_point])
                    rotation = max(min(prev_rotation, np.deg2rad(45)), np.deg2rad(-45))
                    dont_smooth_rotation = True
                    travelling_over_air = True
                elif travelling_over_air:
                    continue
                else:
                    continue
            else:
                if travelling_over_air:
                    if SAFE_TRAVEL_TRANSITIONS:  # not in the notebook
                        reentry = new_position.copy()
                        new_position[2] = max(highest_printed_point, new_position[2])
                    else:
                        new_position[2] = highest_printed_point
                    rotation = max(min(rotation, np.deg2rad(45)), np.deg2rad(-45))
                    dont_smooth_rotation = True
                travelling_over_air = False

            move_command = command
            if reentry is not None and command == "G01" and extrusion is not None and extrusion != RETRACTION_LENGTH and extrusion != -RETRACTION_LENGTH:
                move_command = "G00"; extrusion = None; lost_vertices.append(position)  # not in the notebook
            extrusion_multiplier = 1
            if extrusion is not None and extrusion != RETRACTION_LENGTH and extrusion != -RETRACTION_LENGTH:
                if SMOOTH_EXTRUSION_MULTIPLIER:  # not in the notebook
                    vi, bc = last_bary["v"]
                    extrusion_multiplier = extrusion_multiplier * np.sum(vertex_volume_ratio[vi] * bc)
                else:
                    extrusion_multiplier = extrusion_multiplier * z_squish_scales[containing_cell_index]
                if EXTRUSION_MULTIPLIER_RANGE is not None:  # not in the notebook
                    extrusion_multiplier = min(max(extrusion_multiplier, EXTRUSION_MULTIPLIER_RANGE[0]), EXTRUSION_MULTIPLIER_RANGE[1])
                extrusion = extrusion * min(extrusion_multiplier, MAX_EXTRUSION_MULTIPLIER)
            elif extrusion == -RETRACTION_LENGTH:
                travelling = True
            elif extrusion == RETRACTION_LENGTH:
                travelling = False
            if prev_rotation is not None and not dont_smooth_rotation:
                rotation = ROTATION_AVERAGING_ALPHA * rotation + (1 - ROTATION_AVERAGING_ALPHA) * prev_rotation

            # --- not in the notebook (off by default): SPLIT_RETRACTIONS puts the +/-RETRACTION_LENGTH E value
            # on its own zero-motion line instead of extruding it during the 1 mm travel lift/plunge
            split = SPLIT_RETRACTIONS and extrusion is not None and (extrusion == RETRACTION_LENGTH or extrusion == -RETRACTION_LENGTH)
            motion_extrusion = None if split else extrusion
            if split and extrusion < 0:  # retract in place first, then lift
                here = prev_new_position if prev_new_position is not None else new_position
                new_gcode_points.append({
                    "position": here.copy(), "original_position": position,
                    "rotation": prev_rotation if prev_new_position is not None else rotation,
                    "command": "G01", "extrusion": extrusion, "inv_time_feed": None,
                    "extrusion_multiplier": extrusion_multiplier, "feed": gcode_point["feed"],
                    "travelling": prev_travelling, "e_only": True})

            if prev_rotation is not None and prev_new_position is not None and np.abs(rotation - prev_rotation) > ROTATION_MAX_DELTA:
                delta_rotation = rotation - prev_rotation
                num_interpolations = int(np.abs(delta_rotation) / ROTATION_MAX_DELTA) + 1
                delta_pos = new_position - prev_new_position
                for i in range(num_interpolations):
                    new_gcode_points.append({
                        "position": prev_new_position + (delta_pos * ((i + 1) / num_interpolations)),
                        "original_position": position,
                        "rotation": prev_rotation + (delta_rotation * ((i + 1) / num_interpolations)),
                        "command": move_command if SAFE_TRAVEL_TRANSITIONS else prev_command,
                        "extrusion": motion_extrusion / num_interpolations if motion_extrusion is not None else None,
                        "inv_time_feed": inv_time_feed * num_interpolations if inv_time_feed is not None else None,
                        "extrusion_multiplier": extrusion_multiplier,
                        "feed": gcode_point["feed"],
                        "travelling": prev_travelling
                    })
            else:
                new_gcode_points.append({
                    "position": new_position,
                    "original_position": position,
                    "rotation": rotation,
                    "command": move_command,
                    "extrusion": motion_extrusion,
                    "inv_time_feed": inv_time_feed,
                    "extrusion_multiplier": extrusion_multiplier,
                    "feed": gcode_point["feed"],
                    "travelling": travelling
                })

            if reentry is not None and reentry[2] != new_position[2]:  # not in the notebook: lower as a travel
                new_gcode_points.append({
                    "position": reentry.copy(), "original_position": position, "rotation": rotation,
                    "command": "G00", "extrusion": None, "inv_time_feed": None,
                    "extrusion_multiplier": extrusion_multiplier, "feed": gcode_point["feed"],
                    "travelling": travelling})
                new_position = reentry

            if split and extrusion > 0:  # plunge first, then unretract in place
                if new_gcode_points[-1]["travelling"] != travelling:  # interpolated motion stayed hopped
                    new_gcode_points.append({
                        "position": new_position.copy(), "original_position": position, "rotation": rotation,
                        "command": command, "extrusion": None, "inv_time_feed": None,
                        "extrusion_multiplier": extrusion_multiplier, "feed": gcode_point["feed"],
                        "travelling": travelling})
                new_gcode_points.append({
                    "position": new_position.copy(), "original_position": position, "rotation": rotation,
                    "command": "G01", "extrusion": extrusion, "inv_time_feed": None,
                    "extrusion_multiplier": extrusion_multiplier, "feed": gcode_point["feed"],
                    "travelling": travelling, "e_only": True})

            prev_rotation = rotation
            prev_new_position = new_position.copy()
            prev_travelling = travelling
            prev_command = move_command if SAFE_TRAVEL_TRANSITIONS else command

            if command == "G01" and extrusion is not None and extrusion > 0 and (highest_printed_point != 0 or new_position[2] < 1):
                highest_printed_point = max(highest_printed_point, new_position[2])

    stats = {
        "lost_vertices": len(lost_vertices),
        "gcode_points": len(gcode_points),
        "no_containing_cell": len(no_cell_positions),
        "bad_barycentric_sum": len(bad_barycentric_positions),
    }
    return new_gcode_points, stats


def write_gcode(new_gcode_points, out_path, NOZZLE_OFFSET=42):
    """Cell 18."""
    prev_r = 0
    prev_theta = 0
    prev_z = 20

    theta_accum = 0

    with open(out_path, 'w') as fh:
        fh.write("G94 ; mm/min feed  \n")
        fh.write("G28 ; home \n")
        fh.write("M83 ; relative extrusion \n")
        fh.write("G1 E10 ; prime extruder \n")
        fh.write("G94 ; mm/min feed \n")
        fh.write("G90 ; absolute positioning \n")
        fh.write(f"G0 C{prev_theta} X{prev_r} Z{prev_z} B0 ; go to start \n")
        fh.write("G93 ; inverse time feed \n")

        for i, point in enumerate(new_gcode_points):
            position = point["position"]
            rotation = point["rotation"]

            if np.all(np.isnan(position)):
                continue

            if position[2] < 0:
                continue

            z_hop = 0
            if point["travelling"]:
                z_hop = 1

            r = np.linalg.norm(position[:2])
            theta = np.arctan2(position[1], position[0])
            z = position[2]

            r += -np.sin(rotation) * (NOZZLE_OFFSET + z_hop)
            z += (np.cos(rotation) - 1) * (NOZZLE_OFFSET + z_hop) + z_hop

            delta_theta = theta - prev_theta
            if delta_theta > np.pi:
                delta_theta -= 2 * np.pi
            if delta_theta < -np.pi:
                delta_theta += 2 * np.pi

            theta_accum += delta_theta

            string = f"{point['command']} C{np.rad2deg(theta_accum):.5f} X{r:.5f} Z{z:.5f} B{np.rad2deg(rotation):.5f}"

            if point["extrusion"] is not None:
                string += f" E{point['extrusion']:.4f}"

            no_feed_value = False
            if point.get("e_only"):  # SPLIT_RETRACTIONS: zero-motion retract/unretract at the planar feed
                string += f" F{point['feed']:g}"
                fh.write(f"G94\n")
                no_feed_value = True
            elif point["inv_time_feed"] is not None:
                string += f" F{(point['inv_time_feed']):.4f}"
            else:
                string += f" F20000"
                fh.write(f"G94\n")
                no_feed_value = True

            fh.write(string + "\n")

            if no_feed_value:
                fh.write(f"G93\n")

            prev_r = r
            prev_theta = theta
            prev_z = z
