"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL

Mathematical calculations for lattice element positioning and transformations.
"""

import numpy as np


def transform(x, y, rotation_deg, dx):
    """
    Transform coordinates based on angle and displacement.
    """
    rotation_rad = np.radians(rotation_deg)
    x_new = x + dx * np.cos(rotation_rad)
    y_new = y + dx * np.sin(rotation_rad)
    return x_new, y_new


def rotate_corners(x: float, y: float, rotation_deg: float, ds: float = 1.0, width: float = 0.1) -> np.ndarray:
    """
    Generates rectangle's corners after applying rotation matrix.
    This is utilized to properly visualize a rotated lattice element in Plotly.

    :param x: starting x-coordinate before the rotation
    :param y: starting y-coordinate before the rotation
    :param rotation_deg: Rotation angle in degrees, counterclockwise.
    :param ds: Length of the rectangle along the local X-axis (default is 1.0).
    :param width: Half of the rectangle's height (default is 0.1).
    :return: A NumPy array of shape (5, 2) with rotated (x, y) corner coordinates, closed for polygon plotting.
    """
    rotation_rad = np.radians(rotation_deg)

    corners = np.array([
        [0, -width],
        [ds, -width],
        [ds, width],
        [0, width],
        [0, -width]  # close polygon
    ])

    R = np.array([
        [np.cos(rotation_rad), -np.sin(rotation_rad)],
        [np.sin(rotation_rad),  np.cos(rotation_rad)],
    ])

    rotated = corners @ R.T + [x, y]
    return rotated
