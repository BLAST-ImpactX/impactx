"""
The 3D survey: the reference orbit and the cameras that follow it, see Survey3DGeometry
for the geometry of the lattice

This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""

import math

#: element kinds that are not drawn, as in plot_survey over s
SKIP_KINDS = ("Drift", "ChrDrift", "ExactDrift", "Empty", "Marker", "Source")
#: elements drawn as beam pipes in lattices of opaque elements
PIPE_KINDS = ("Drift", "ChrDrift", "ExactDrift")

#: largest turn of the reference trajectory per cell of the geometry [rad]
MAX_CELL_ANGLE = math.radians(2.5)

#: longitudinal wires of the envelope cage
_CAGE_WIRES = 8

#: vertices around a ring of the envelope cage
_RING_VERTICES = 16

#: names in the Blueprint mesh of the geometry: they differ from the ones of the beam
SURVEY_TOPOLOGY = "survey"
SURVEY_FIELD = "survey_kind"
WIRES_TOPOLOGY = "survey_wires"
WIRES_FIELD = "survey_wire_kind"

#: values of WIRES_FIELD
WIRE_ORBIT = 0.0
WIRE_ENVELOPE = 1.0


#: color of the envelope tube
ENVELOPE_COLOR = (0.95, 0.55, 0.15)


def frames(df):
    """Position and local axes of the reference particle per survey row

    The local axes are the direction of motion ``t``, the horizontal axis ``n`` of
    the particle coordinate x, and the vertical axis ``v`` of y. The frame has no
    roll: ``n`` stays horizontal (in the global x-z plane).

    Returns
    -------
    tuple of numpy.ndarray
        p, t, n, v, each of shape (rows, 3)

    Raises
    ------
    ValueError
        If the reference particle moves vertically or does not move.
    """
    import numpy as np

    p = df[["x", "y", "z"]].to_numpy(dtype=np.float64)
    mom = df[["px", "py", "pz"]].to_numpy(dtype=np.float64)
    norm = np.linalg.norm(mom, axis=1)
    if np.any(norm == 0.0):
        raise ValueError("survey 3D: the reference particle does not move.")
    t = mom / norm[:, None]
    up = np.array([0.0, 1.0, 0.0])
    n = np.cross(up, t)
    n_norm = np.linalg.norm(n, axis=1)
    if np.any(n_norm < 1.0e-12):
        raise ValueError("survey 3D: the reference particle moves vertically.")
    n /= n_norm[:, None]
    v = np.cross(t, n)
    return p, t, n, v


def overview_axis(df):
    """Direction in the x-z plane that points to the right in an overview of the lattice

    This is the long axis of the orbit, along the net travel of the reference particle.
    Without net travel along it, as in a ring, it follows the direction of motion at the
    start of the lattice. Without a long axis, as in a ring of equal cells, it is that
    direction of motion. The result depends on the orbit only, so that an overview stays
    in place while the beam moves.

    Parameters
    ----------
    df : pandas.DataFrame
        The survey of the lattice.

    Returns
    -------
    numpy.ndarray
        Unit vector (x, z).
    """
    import numpy as np

    # relative: well above the roundoff of a survey, in single precision too
    tolerance = 1.0e-3

    p, t, _, _ = frames(df)
    start = t[0, [0, 2]] / np.linalg.norm(t[0, [0, 2]])
    xz = p[:, [0, 2]] - p[:, [0, 2]].mean(axis=0)
    if np.allclose(xz, 0.0):
        return start

    singular, axes = np.linalg.svd(xz, full_matrices=False)[1:]
    axis = axes[0] if singular[0] - singular[1] > tolerance * singular[0] else start

    size = np.linalg.norm(xz, axis=1).max()
    for forward in (axis @ (xz[-1] - xz[0]) / size, axis @ start, axis @ xz[0] / size):
        if abs(forward) > tolerance:
            return axis if forward > 0.0 else -axis
    return axis


def _closed(p, t, s):
    """Does the orbit close on itself, as in a ring?"""
    import numpy as np

    length = s[-1] - s[0]
    return bool(
        length > 0.0
        and np.linalg.norm(p[-1] - p[0]) <= 1.0e-3 * length
        and np.dot(t[-1], t[0]) >= np.cos(1.0e-3)
    )


def orbit_position(df, s0):
    """Position on the reference orbit at a position s

    Interpolated linearly between the rows of the survey. If the orbit closes on
    itself, as in a ring, it wraps around. Otherwise, the orbit continues straight
    before and after the lattice.

    Parameters
    ----------
    df : pandas.DataFrame
        A survey of the lattice.
    s0 : float
        Position along the orbit [m].

    Returns
    -------
    numpy.ndarray
        The position x, y, z [m].
    """
    import numpy as np

    p, t, _, _ = frames(df)
    s = df["s"].to_numpy(dtype=np.float64)
    if _closed(p, t, s):
        s0 = s[0] + (s0 - s[0]) % (s[-1] - s[0])
    if s0 <= s[0]:
        return p[0] + (s0 - s[0]) * t[0]
    if s0 >= s[-1]:
        return p[-1] + (s0 - s[-1]) * t[-1]
    i = int(np.searchsorted(s, s0, side="right"))
    # the rows around s0, s[i - 1] <= s0 < s[i]: zero-length elements share an s
    h = s[i] - s[i - 1]
    w = (s0 - s[i - 1]) / h if h > 0.0 else 0.0
    return (1.0 - w) * p[i - 1] + w * p[i]


def smoothed_axes(df, s0, sigma):
    """Axes of the reference orbit, averaged with a Gaussian weight around a position

    A camera that follows the beam along these axes turns steadily through a lattice
    of bends, instead of turning in each bend and not between them. If the orbit
    closes on itself, as in a ring, the average wraps around. Otherwise, the orbit
    continues straight before and after the lattice.

    Parameters
    ----------
    df : pandas.DataFrame
        A survey of the lattice.
    s0 : float
        Center of the average [m].
    sigma : float
        Width of the Gaussian weight in s [m], 0 for the axes at s0.

    Returns
    -------
    tuple of numpy.ndarray
        The direction of motion t, the horizontal axis n and the vertical axis v.
    """
    import numpy as np

    p, t, _, _ = frames(df)
    s = df["s"].to_numpy(dtype=np.float64)
    ds = s - s0
    length = s[-1] - s[0]
    closed = _closed(p, t, s)
    if closed:
        ds = (ds + 0.5 * length) % length - 0.5 * length
    local = t[np.argmin(np.abs(ds))]
    if sigma > 0.0 and len(s) > 1:
        # each row weighted by the length of orbit that it stands for (trapezoids)
        h = np.abs(np.diff(s))
        dl = np.zeros_like(s)
        dl[:-1] += 0.5 * h
        dl[1:] += 0.5 * h
        weight = np.exp(-0.5 * (ds / sigma) ** 2) * dl
        average = weight @ t
        if not closed:
            # the straight orbit before and after the lattice: Gaussian tails
            tail = sigma * math.sqrt(0.5 * math.pi)
            before = tail * math.erfc((s0 - s[0]) / (math.sqrt(2.0) * sigma))
            after = tail * math.erfc((s[-1] - s0) / (math.sqrt(2.0) * sigma))
            average = average + before * t[0] + after * t[-1]
        horizontal = np.linalg.norm(average[[0, 2]])
        # e.g., a window over a whole ring averages out
        if horizontal > 1.0e-6 * sigma:
            local = average / np.linalg.norm(average)
    up = np.array([0.0, 1.0, 0.0])
    n = np.cross(up, local)
    n /= np.linalg.norm(n)
    return local, n, np.cross(local, n)


def refined_lattice(lattice, df, max_length=None, pipes=False):
    """Copies of the elements with more slices, for a smooth geometry

    Each slice turns the reference trajectory by at most MAX_CELL_ANGLE and, for
    elements that are drawn, is at most max_length long. The number of slices of an
    element is a multiple of its own, so its slice boundaries stay slice boundaries.

    Parameters
    ----------
    lattice : KnownElementsList
    df : pandas.DataFrame
        The survey of the lattice, one row per slice.
    max_length : float, optional
        Longest slice of a drawn element [m].
    pipes : bool
        Drifts are drawn, as beam pipes.

    Returns
    -------
    KnownElementsList, or the lattice itself if no element needs more slices
    """
    import numpy as np

    from ..extensions.KnownElementsList import _element_kind
    from ..impactx_pybind import elements

    _, t, _, _ = frames(df)
    element = df["element"].to_numpy()
    # turn of the direction of motion in each slice: from the row before to its row
    cos_turn = np.clip(np.sum(t[1:] * t[:-1], axis=1), -1.0, 1.0)
    turn = np.arccos(cos_turn)

    # the largest turn per element
    largest_turn = np.zeros(len(lattice))
    np.maximum.at(largest_turn, element[1:], turn)

    factors = {}
    for index in np.unique(element[1:]):
        el = lattice[int(index)]
        # each read of a property is a call into C++
        ds = el.ds
        if ds == 0.0:
            continue
        largest = float(largest_turn[index])
        factor = math.ceil(largest / MAX_CELL_ANGLE)
        kind = _element_kind(el)
        drawn = kind not in SKIP_KINDS or (pipes and kind in PIPE_KINDS)
        if max_length is not None and drawn:
            factor = max(factor, math.ceil(abs(ds) / el.nslice / max_length))
        if factor > 1:
            factors[int(index)] = factor
    if not factors:
        return lattice

    refined = elements.KnownElementsList()
    for index, el in enumerate(lattice):
        if index in factors:
            el = el.copy()
            el.nslice = el.nslice * factors[index]
        refined.append(el)
    return refined


def look_from(points, direction, up, aspect, fov_deg, margin=1.05):
    """Camera position that looks at all points from a direction, in perspective

    Parameters
    ----------
    points : numpy.ndarray
        Points to fit into the view, of shape (num, 3).
    direction : numpy.ndarray
        From the look-at point to the camera.
    up : numpy.ndarray
        Up direction of the image.
    aspect : float
        Image width / height.
    fov_deg : float
        Vertical field of view [deg].

    Returns
    -------
    dict
        Ascent camera parameters.
    """
    import numpy as np

    w = np.asarray(direction, dtype=np.float64)
    w = w / np.linalg.norm(w)
    right = np.cross(up, w)
    right /= np.linalg.norm(right)
    up_c = np.cross(w, right)

    center = 0.5 * (points.min(axis=0) + points.max(axis=0))
    q = points - center
    tan_v = math.tan(math.radians(fov_deg) / 2.0)
    tan_h = aspect * tan_v
    depth = q @ w
    distance = np.max(
        np.maximum(
            depth + np.abs(q @ right) * margin / tan_h,
            depth + np.abs(q @ up_c) * margin / tan_v,
        )
    )
    return {
        "look_at": center.tolist(),
        "position": (center + distance * w).tolist(),
        "up": up_c.tolist(),
        "fov": fov_deg,
    }
