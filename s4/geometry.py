"""
Small numpy geometry shared by the deformation, the mapper and the G-code checks.

Machine frame: the bed rotates by C (deg) about z; X is the radial position (mm) and Z the height of the B pivot;
B (deg) tilts the nozzle in the radial-vertical plane, about the tangential axis; L is the pivot-to-tip distance
(NOZZLE_OFFSET). fast_map.write_gcode turns a tip at radius r, height z into
    X = r - sin(B) * L,    Z = z + (cos(B) - 1) * L
(plus the travel z-hop), so nozzle_tip inverts that.
"""
import numpy as np


def tangential_vectors(cell_centers):
    """Unit vector along +z x (x, y, 0) per cell centre: the axis the B tilt rotates about there ([1, 0, 0] on the
    z axis, where it is undefined)."""
    cxy = cell_centers[:, :2]
    c3 = np.hstack([cxy, np.zeros((cxy.shape[0], 1))])
    t = np.cross(np.array([0, 0, 1]), c3)
    with np.errstate(invalid="ignore"):
        t /= np.linalg.norm(t, axis=1)[:, None]
    t[np.isnan(t).any(axis=1)] = [1, 0, 0]
    return t


def rotation_matrices(cell_centers, rotation_field):
    """3x3 rotation per cell: rotation_field (rad) about the cell's tangential axis."""
    from scipy.spatial.transform import Rotation as R
    return R.from_rotvec(rotation_field[:, None] * tangential_vectors(cell_centers)).as_matrix()


def nozzle_tip(C, X, Z, B, L):
    """Nozzle-tip position (x, y, z) in the bed frame for machine axes C, X, Z, B (deg, mm, mm, deg) and nozzle
    offset L, as a tuple of three. Scalars or arrays (elementwise)."""
    b = np.radians(B)
    th = np.radians(C)
    r = X + np.sin(b) * L
    z = Z - (np.cos(b) - 1) * L
    return r * np.cos(th), r * np.sin(th), z
