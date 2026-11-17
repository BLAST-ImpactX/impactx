"""
Meshes of the 3D survey: elements swept along the reference orbit, as hexahedra or
quadrilaterals, and their lighting helpers

This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""


def _shade(rgb, factor, tint=(0.0, 0.0, 0.0)):
    """A color times a brightness factor, plus a tint"""
    return tuple(min(factor * c + d, 1.0) for c, d in zip(rgb, tint))


def _smoothstep(edge0, edge1, x):
    """0 below edge0, 1 above edge1, smooth in between"""
    import numpy as np

    t = np.clip((x - edge0) / (edge1 - edge0), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _row_dots(a, b):
    """Dot products of the rows of two arrays of shape (num, 3)"""
    import numpy as np

    # faster than np.sum(a * b, axis=1)
    return np.einsum("ij,ij->i", a, b)


def _row_norms(vectors):
    """Lengths of the rows of an array of shape (num, 3)"""
    import numpy as np

    # faster than np.linalg.norm(vectors, axis=1)
    return np.sqrt(_row_dots(vectors, vectors))


def _normalized(vectors):
    """Unit vectors, of shape (num, 3)"""
    return vectors / _row_norms(vectors)[:, None]


def _unchanged(old, new, tolerance):
    """Are two arrays of the same shape, and equal within an absolute tolerance?"""
    import numpy as np

    return old.shape == new.shape and (
        old.size == 0 or float(np.max(np.abs(old - new))) <= tolerance
    )


def _unique_positions(s, rows):
    """Rows with distinct s: the last one at each position"""
    import numpy as np

    rows = np.asarray(rows)
    keep = np.ones(len(rows), dtype=bool)
    keep[:-1] = s[rows[1:]] != s[rows[:-1]]
    return rows[keep]


def _add_hexahedra(mesh, name, points, cells, shape="hex"):
    """Add an unstructured topology of hexahedra (or of another shape) and its
    coordinate set to a mesh, as numpy arrays"""
    import numpy as np

    coordset = f"{name}_coords"
    mesh["coordsets"][coordset] = {
        "type": "explicit",
        "values": {
            axis: np.ascontiguousarray(points[:, k], dtype=np.float64)
            for k, axis in enumerate("xyz")
        },
    }
    mesh["topologies"][name] = {
        "type": "unstructured",
        "coordset": coordset,
        "elements": {
            "shape": shape,
            "connectivity": np.ascontiguousarray(cells.ravel(), dtype=np.int64),
        },
    }


def _surface_patches(
    p, t, n, v, half_width, half_height, radius, bore, outside="box", lower_half=False
):
    """The surfaces of a box with a bore, or of its lower half, cut open at y = 0

    Flat faces, i.e., the sides of the box and of a square bore, the cut and the end
    caps, are separate patches with their own vertices, so that their edges stay
    sharp; a round bore is one patch with radial normals, smooth.

    Parameters
    ----------
    p, t, n, v : numpy.ndarray
        Positions and axes of the frames along the element, (frames, 3).
    radius : float
        Radius of a round bore, half size of a square bore, 0 for a solid box.
    bore : str
        ``"round"`` or ``"square"``.
    outside : str
        ``"box"``, or ``"round"`` for a pipe of the radius half_width, one smooth patch.
    lower_half : bool
        Only the visible surfaces of the lower half (y <= 0): cut open through the axis.

    Returns
    -------
    The vertices (num, 3), the quadrilaterals (num, 4) and the normals per vertex.
    """
    import numpy as np

    hw, hh = half_width, half_height
    if lower_half:
        # rays from the axis over the lower half: uniform, and to the corners of box
        # and bore
        corners = np.mod(np.arctan2([-hh, -hh], [-hw, hw]), 2.0 * np.pi)
        uniform = np.pi + np.pi * np.arange(33) / 32
        phi = np.unique(
            np.round(
                np.concatenate((uniform, corners, [1.25 * np.pi, 1.75 * np.pi])), 12
            )
        )
    else:
        # rays from the axis all around: uniform (including the corners of a square
        # bore at 45 degrees), and to the corners of the box; closed, the first ray
        # again at the end
        corners = np.arctan2([-hh, -hh, hh, hh], [-hw, hw, hw, -hw])
        uniform = 2.0 * np.pi * np.arange(64) / 64
        phi = np.unique(
            np.round(np.mod(np.concatenate((uniform, corners)), 2.0 * np.pi), 12)
        )
        phi = np.append(phi, phi[0] + 2.0 * np.pi)
    c, s = np.cos(phi), np.sin(phi)
    c[np.abs(c) < 1.0e-12] = 0.0
    s[np.abs(s) < 1.0e-12] = 0.0
    with np.errstate(divide="ignore"):
        reach = np.minimum(hw / np.abs(c), hh / np.abs(s))
        reach_bore = np.minimum(1.0 / np.abs(c), 1.0 / np.abs(s))
    if outside == "round":
        outer = np.stack((hw * c, hw * s), axis=1)
    else:
        outer = np.stack((reach * c, reach * s), axis=1)
    scale = radius * (reach_bore if bore == "square" else 1.0)
    inner = np.stack((scale * c, scale * s), axis=1)

    parts = []

    def strip(section, normals2d):
        """A patch along the frames, through the points of a cross section"""
        m = len(section)
        points = (
            p[:, None, :]
            + section[None, :, 0, None] * n[:, None, :]
            + section[None, :, 1, None] * v[:, None, :]
        )
        normals = (
            normals2d[None, :, 0, None] * n[:, None, :]
            + normals2d[None, :, 1, None] * v[:, None, :]
        )
        base = m * np.arange(len(p) - 1)[:, None] + np.arange(m - 1)[None, :]
        quads = np.stack((base, base + 1, base + m + 1, base + m), axis=-1)
        parts.append(
            (points.reshape(-1, 3), quads.reshape(-1, 4), normals.reshape(-1, 3))
        )

    def flat_sides(section, half_x, sign):
        """Patches of the sides of a rectangle, with their normals (outward: sign 1)"""
        middle = 0.5 * (section[:-1] + section[1:])
        on_x = np.isclose(np.abs(middle[:, 0]), half_x)
        normal = np.where(
            on_x[:, None],
            np.stack((np.sign(middle[:, 0]), np.zeros(len(middle))), axis=1),
            np.stack((np.zeros(len(middle)), np.sign(middle[:, 1])), axis=1),
        )
        start = 0
        for i in range(1, len(middle) + 1):
            if i == len(middle) or not np.array_equal(normal[i], normal[start]):
                side = section[start : i + 1]
                strip(side, np.tile(sign * normal[start], (len(side), 1)))
                start = i

    if outside == "round":
        strip(outer, np.stack((c, s), axis=1))
    else:
        flat_sides(outer, hw, 1.0)
    if radius > 0.0:
        if bore == "square":
            flat_sides(inner, radius, -1.0)
        else:
            strip(inner, -np.stack((c, s), axis=1))
    if lower_half:
        # the cut, facing up
        up = np.tile([0.0, 1.0], (2, 1))
        strip(np.array([outer[0], inner[0]]), up)
        strip(np.array([inner[-1], outer[-1]]), up)
    # the end caps
    m = len(phi)
    for j, normal in ((0, -t[0]), (len(p) - 1, t[-1])):
        ring = np.concatenate((outer, inner))
        points = p[j] + ring[:, 0, None] * n[j] + ring[:, 1, None] * v[j]
        i = np.arange(m - 1)
        quads = np.stack((i, i + 1, m + i + 1, m + i), axis=-1)
        parts.append((points, quads, np.tile(normal, (len(points), 1))))

    points, quads, normals = [], [], []
    offset = 0
    for part_points, part_quads, part_normals in parts:
        points.append(part_points)
        quads.append(part_quads + offset)
        normals.append(part_normals)
        offset += len(part_points)
    return np.concatenate(points), np.concatenate(quads), np.concatenate(normals)


def _sweep(p, n, v, section, quads):
    """Sweep a cross section along frames into hexahedra

    Parameters
    ----------
    p, n, v : numpy.ndarray
        Positions and transverse axes of the frames, (frames, 3).
    section : numpy.ndarray
        Vertices of the cross section per frame, (frames, vertices, 2) in (n, v).
    quads : numpy.ndarray
        Convex quadrilaterals of the cross section, counter-clockwise around the
        direction of motion, as indices into its vertices, (quads, 4).

    Returns the vertices (frames * vertices, 3) and the hexahedra
    ((frames - 1) * quads, 8).
    """
    import numpy as np

    points = (
        p[:, None, :]
        + section[:, :, 0, None] * n[:, None, :]
        + section[:, :, 1, None] * v[:, None, :]
    ).reshape(-1, 3)
    per_frame = section.shape[1]
    base = per_frame * np.arange(len(p) - 1)
    cells = np.concatenate(
        (
            base[:, None, None] + quads[None, :, :],
            base[:, None, None] + per_frame + quads[None, :, :],
        ),
        axis=2,
    ).reshape(-1, 8)
    return points, cells


def _box_sweep_varying(p, n, v, half_width, half_height):
    """A box swept along frames, with a half width and height per frame"""
    import numpy as np

    corners = np.array(((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)))
    section = (
        corners[None, :, :] * np.stack((half_width, half_height), axis=1)[:, None, :]
    )
    return _sweep(p, n, v, section, np.array([[0, 1, 2, 3]]))


def _shell_sweep(
    p,
    n,
    v,
    half_width,
    half_height,
    radius,
    bore="round",
    liner=0.0,
    lower_half=False,
    seam=0.0,
):
    """A box with a bore swept along frames

    The cross section is split along rays from the center into convex
    quadrilaterals, between the outer rectangle and the bore. The rays include the
    corners of the rectangle, so that each quadrilateral has its outer edge on one
    side of the rectangle and of a square bore. With a liner, a thin layer of
    quadrilaterals lines the bore.

    Parameters
    ----------
    radius : float
        Radius of a round bore, half size of a square bore.
    bore : str
        ``"round"`` or ``"square"``.
    liner : float
        Thickness of the lining of the bore, 0 for none.
    lower_half : bool
        Only the lower half (y <= 0), cut open through the axis.
    seam : float
        With a liner: thickness of the layer that joins the wall and the lining.

    Returns
    -------
    The vertices, the hexahedra, and per vertex if it lines the bore and its normal:
    radial, towards the axis for the lining.
    """
    import numpy as np

    corners = np.arctan2(
        [-half_height, -half_height, half_height, half_height],
        [-half_width, half_width, half_width, -half_width],
    )
    # rays: uniform (including the corners of a square bore at 45 degrees) and to the
    # corners of the box
    uniform = 2.0 * np.pi * np.arange(32) / 32
    phi = np.unique(
        np.round(np.mod(np.concatenate((uniform, corners)), 2.0 * np.pi), 12)
    )
    if lower_half:
        # from the cut at -x along the bottom to the cut at +x
        phi = np.concatenate(([np.pi], phi[np.sin(phi) < -1.0e-12], [2.0 * np.pi]))
    c, s = np.cos(phi), np.sin(phi)
    if bore == "square":
        # the square of half size 1 along the rays
        with np.errstate(divide="ignore"):
            reach_unit = np.minimum(1.0 / np.abs(c), 1.0 / np.abs(s))
        unit = np.stack((reach_unit * c, reach_unit * s), axis=1)
    else:
        unit = np.stack((c, s), axis=1)
    # rings from the outside in: box, wall to the seam, seam to the lining, lining,
    # bore. Wall and lining have their own vertices, the seam joins them: their shading
    # is interpolated over the vertices of each cell.
    if liner > 0.0:
        rings_radii = [radius + liner + seam, radius + liner, radius]
        ring_lined = [False, False, True, True]
    else:
        rings_radii = [radius]
        ring_lined = [False, False]
    # the outer rectangle along the same rays
    with np.errstate(divide="ignore"):
        reach = np.minimum(half_width / np.abs(c), half_height / np.abs(s))
    outer = np.stack((reach * c, reach * s), axis=1)

    # rings of vertices from the outside in: box, lining, bore
    rings = [outer] + [r * unit for r in rings_radii]
    m = len(phi)
    section = np.broadcast_to(np.concatenate(rings), (len(p), len(rings) * m, 2))
    # counter-clockwise: along the outer ring, then back along the inner ring; a cut
    # open half does not close around
    sectors = m - 1 if lower_half else m
    quads = np.array(
        [
            [k * m + i, k * m + (i + 1) % m, (k + 1) * m + (i + 1) % m, (k + 1) * m + i]
            for k in range(len(rings) - 1)
            for i in range(sectors)
        ]
    )
    points, cells = _sweep(p, n, v, section, quads)

    # per vertex: lining or wall, and the normal: radial, outward for the walls and
    # towards the axis for the lining. In translucent elements, the walls are seen
    # through: the face normals of the box would show the borders between the faces,
    # along the diagonals.
    radial = unit / np.linalg.norm(unit, axis=1)[:, None]
    section_lined = np.repeat(ring_lined, m)
    section_normals = np.concatenate(
        [-radial if lined else radial for lined in ring_lined]
    )
    lined = np.tile(section_lined, len(p))
    normals = (
        section_normals[None, :, 0, None] * n[:, None, :]
        + section_normals[None, :, 1, None] * v[:, None, :]
    ).reshape(-1, 3)
    return points, cells, lined, normals
