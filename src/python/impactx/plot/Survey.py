"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def survey(self, ref, *, slices=True, covariance=None, element=None):
    """Survey of the lattice: the reference trajectory in the global frame.

    The reference particle can be at any slice boundary of the lattice, e.g., the reference
    particle of a beam during or after tracking. A copy of it is pushed backward to the start
    of the lattice, then forward through the elements, slice by slice, as in tracking. The
    reference particle itself is not changed.

    With the covariance matrix of the beam at the position of the reference particle, the
    survey also contains the rms beam sizes along the lattice, from envelope tracking of the
    covariance matrix: backward to the start of the lattice and forward from there. This
    needs a linear transport map for each element.

    Parameters
    ----------
    self : KnownElementsList
        The lattice.
    ref : RefPart
        The reference particle at a slice boundary of the lattice, e.g., at its start. In a
        repeated lattice, its position is taken within the current period.
    slices : bool
        One row after each element slice if true, otherwise one row after each element.
    covariance : Map6x6, optional
        Covariance matrix of the beam at the position of the reference particle, e.g., from
        ``ImpactXParticleContainer.covariance_matrix()``.
    element : optional
        The element of this lattice, or its index, that the reference particle is in or has
        just passed, e.g., ``sim.tracking_element`` in a tracking hook. It is needed if the
        reference particle is inside the lattice at the position of zero-length elements that
        change the reference particle or the beam. At the start of the lattice, the reference
        particle is before such elements, at the end after them. An element that is in the
        lattice more than once is the one at the position of the reference particle.

    Returns
    -------
    pandas.DataFrame
        One row for the start and one after each slice or element, with the columns:

        - ``element``: index of the element in the lattice, ``-1`` for the start
        - ``name``, ``kind``: name and kind of the element
        - ``s``: path length of the reference trajectory [m]
        - ``x``, ``y``, ``z``: position in the global frame [m]
        - ``theta``: horizontal angle of the direction of motion to the :math:`z` axis,
          towards :math:`+x` [rad]
        - ``phi``: vertical angle of the direction of motion to the :math:`x`-:math:`z`
          plane [rad]
        - ``px``, ``py``, ``pz``, ``pt``: momenta and energy of the reference particle
        - ``kin_energy_MeV``: kinetic energy of the reference particle [MeV]

        and with a covariance matrix: ``sigma_x``, ``sigma_px``, ``sigma_y``, ``sigma_py``,
        ``sigma_t``, ``sigma_pt``, the rms beam sizes in the units of the phase space
        coordinates.

    Raises
    ------
    ValueError
        If the reference particle is not at a slice boundary of the lattice or of the
        element, if its position among zero-length elements is ambiguous without the
        element, or with an element that is in the lattice more than once and is not its
        index, or if the element is not in the lattice.
    RuntimeError
        With a covariance matrix, if an element has no linear transport map, or the
        reference particle is not at the start and an element cannot be reversed.
    """
    return _survey(self, ref, slices=slices, covariance=covariance, element=element)[0]


def _survey(self, ref, *, slices=True, covariance=None, element=None):
    """survey(), and the index of the element of the reference particle in the lattice

    Returns the survey and the index, None without an element.
    """
    import math

    import pandas as pd

    from .. import Map6x6, push
    from ..extensions.KnownElementsList import _element_kind

    lattice = list(self)
    # properties of the elements, read once: each read is a call into C++
    nslices = [el.nslice for el in lattice]

    # the slices of the lattice, and the path length at their boundaries
    steps = []  # (element index, slice index)
    first_steps = []  # per element: its first slice
    s_grid = [0.0]
    for index, el in enumerate(lattice):
        nslice = nslices[index]
        slice_ds = el.ds / nslice
        first_steps.append(len(steps))
        for i in range(nslice):
            steps.append((index, i))
            s_grid.append(s_grid[-1] + slice_ds)
    length = s_grid[-1]
    tol = 1.0e-9 * max(length, 1.0)

    # position of the reference particle in the current period
    s_local = ref.s
    if length > 0.0 and s_local > length + tol:
        s_local -= math.floor((s_local - tol) / length) * length
    s_start = ref.s - s_local

    def is_identity(el):
        """A zero-length element that changes neither the reference particle nor the beam"""
        probe = ref.copy()
        push(probe, el)
        coords = ("x", "y", "z", "px", "py", "pz", "pt")
        if any(getattr(probe, c) != getattr(ref, c) for c in coords):
            return False
        try:
            R = el.transfer_map(ref)
        except RuntimeError:
            return False
        return all(
            R[i, j] == (1.0 if i == j else 0.0)
            for i in range(1, 7)
            for j in range(1, 7)
        )

    def boundaries(index):
        """The slice boundaries of an element: indices into s_grid"""
        return range(first_steps[index], first_steps[index] + nslices[index] + 1)

    # the element of the reference particle: an element that is in the lattice more
    # than once is the one at the position of the reference particle
    element_index = element
    if element is not None and not isinstance(element, int):
        indices = [i for i, el in enumerate(lattice) if el is element]
        if not indices:
            raise ValueError("survey: the element is not in this lattice.")
        at_ref = [
            i
            for i in indices
            if any(abs(s_grid[k] - s_local) <= tol for k in boundaries(i))
        ]
        element_index = at_ref[-1] if at_ref else indices[0]
        if len(at_ref) > 1:
            # the reference particle is at the end of one and the start of the next:
            # unambiguous only if the elements between them change nothing
            between = [
                lattice[steps[k][0]]
                for k in range(
                    first_steps[at_ref[0]] + nslices[at_ref[0]],
                    first_steps[element_index],
                )
            ]
            if not all(is_identity(el) for el in between):
                raise ValueError(
                    f"survey: the element '{element.name or ''}' is in the lattice more "
                    f"than once at s = {ref.s} m, around zero-length elements that change "
                    "the beam. Pass its index instead."
                )

    # slice boundary of the reference particle: an index into s_grid
    if element_index is None:
        candidates = range(len(s_grid))
    else:
        candidates = boundaries(element_index)
    matches = [k for k in candidates if abs(s_grid[k] - s_local) <= tol]
    if not matches:
        raise ValueError(
            f"survey: the reference particle at s = {ref.s} m is not at a slice "
            "boundary of the lattice (or of the element)."
        )

    # zero-length elements at the position of the reference particle: before or after them
    boundary = matches[-1]
    if len(matches) > 1 and element is None:
        if s_local <= tol:
            boundary = matches[0]  # before the lattice, e.g., before tracking
        elif s_local < length - tol:
            changing = [
                lattice[steps[k][0]]
                for k in range(matches[0], matches[-1])
                if not is_identity(lattice[steps[k][0]])
            ]
            if changing:
                names = ", ".join(
                    f"'{el.name or ''}' ({_element_kind(el)})" for el in changing
                )
                raise ValueError(
                    f"survey: the reference particle at s = {ref.s} m is at the position of "
                    f"zero-length elements that change the beam: {names}. Pass the element "
                    "that it is in or has just passed, e.g., element=sim.tracking_element."
                )

    r = ref.copy()
    cm = Map6x6(covariance) if covariance is not None else None

    def check_linear(el):
        """Elements without a linear transport map raise, before envelope tracking"""
        try:
            el.transfer_map(r)
        except (AttributeError, RuntimeError) as e:
            raise RuntimeError(
                f"survey: envelope tracking through element '{el.name or ''}' "
                f"({_element_kind(el)}): {e}"
            ) from e

    # backward to the start of the lattice, through reversed elements
    reversed_elements = {}
    previous = None
    for k in range(boundary - 1, -1, -1):
        index, i = steps[k]
        el = lattice[index]
        if index not in reversed_elements:
            if not hasattr(el, "reverse"):
                raise RuntimeError(
                    f"survey: element '{el.name or ''}' ({_element_kind(el)}) cannot be "
                    "reversed, to start the survey before the reference particle."
                )
            rev = el.copy()
            rev.reverse()
            if cm is not None:
                check_linear(rev)
            reversed_elements[index] = rev
        if index != previous:
            # the edge of an element is at its entrance in the direction of motion
            r.sedge = s_start + s_grid[k + 1 + (nslices[index] - 1 - i)]
            previous = index
        push(r, reversed_elements[index])
        if cm is not None:
            reversed_elements[index].push(cm, r)

    columns = [
        "element",
        "name",
        "kind",
        "s",
        "x",
        "y",
        "z",
        "theta",
        "phi",
        "px",
        "py",
        "pz",
        "pt",
        "kin_energy_MeV",
    ]
    if cm is not None:
        columns += [f"sigma_{c}" for c in ("x", "px", "y", "py", "t", "pt")]

    def row(index, name, kind):
        # each property of the reference particle is read once
        px, py, pz = r.px, r.py, r.pz
        p = math.sqrt(px**2 + py**2 + pz**2)
        values = (
            index,
            name,
            kind,
            r.s,
            r.x,
            r.y,
            r.z,
            math.atan2(px, pz),
            math.asin(py / p) if p > 0.0 else 0.0,
            px,
            py,
            pz,
            r.pt,
            r.kin_energy_MeV,
        )
        if cm is not None:
            values += tuple(math.sqrt(max(cm[j, j], 0.0)) for j in range(1, 7))
        return values

    # forward through the lattice
    rows = [row(-1, "", "")]
    for index, el in enumerate(lattice):
        name = el.name or ""
        kind = _element_kind(el)
        nslice = nslices[index]
        if cm is not None:
            check_linear(el)
        r.sedge = r.s
        for i in range(nslice):
            push(r, el)
            if cm is not None:
                el.push(cm, r)
            if slices or i == nslice - 1:
                rows.append(row(index, name, kind))

    return pd.DataFrame(rows, columns=columns), element_index


def plot_survey(
    self,
    ref=None,
    ax=None,
    legend=True,
    legend_ncols=5,
    palette="cern-lhc",
    *,
    projection="s",
    beam=None,
    element=None,
    envelope=None,
    magnify="auto",
    camera="overview",
    camera_smoothing="auto",
    image_size=(1280, 720),
    lighting=None,
    bore="round",
    cutaway=False,
    opacity=0.3,
    supersample=2,
    options=None,
):
    """Plot the elements of the KnownElementsList.

    With ``projection="s"``, the elements are drawn over s.
    A positive element strength denotes horizontal focusing (e.g. for quadrupoles) and bending to the right (for dipoles).  In general, this depends on both the sign of the field and the sign of the charge.

    With ``projection="3D"``, Ascent renders the lattice in the global frame of the
    survey (see :py:meth:`survey`): the elements, the reference orbit, and with a
    beam, the beam particles and its envelope.
    The beam particles and the envelope are magnified, to be visible next to the
    elements. The elements are drawn as boxes of a size relative to the lattice.
    This needs ImpactX with Ascent, and is MPI-collective with a beam.

    Parameters
    ----------
    self : KnownElementsList
        The KnownElementsList class in ImpactX
    ref : RefPart
        A reference particle, checked for the charge sign to plot focusing/defocusing strength directions properly.
        For ``projection="3D"`` without a beam, the reference particle at a slice boundary of the lattice.
    ax : matplotlib axes
        A plotting area in matplotlib (called axes there).
    legend: bool
        Plot a legend if true.
    legend_ncols: int
        Number of columns for lattice element types in the legend.
    palette: string
        Color palette.
    projection : str
        ``"s"``: the elements over s. ``"3D"``: rendered with Ascent.
    beam : ImpactXParticleContainer, optional
        For ``projection="3D"``: the particle beam, at a slice boundary of the lattice.
    element : optional
        For ``projection="3D"``: the element of this lattice, or its index, that the
        beam is in or has just passed, see :py:meth:`survey`.
    envelope : float, optional
        For ``projection="3D"``: the half size of the envelope tube in rms beam sizes,
        ``0`` for none. Default: 1 with a beam.
        The envelope needs a linear transport map for each element.
    magnify : str or float or tuple of float
        For ``projection="3D"``: magnification of the beam coordinates and of the
        envelope. ``"auto"``: three rms beam sizes fill most of an element, and the
        bunch is at most a few elements long. A factor, or a pair of factors for the
        transverse and the longitudinal coordinates.
    camera : str or dict
        For ``projection="3D"``: ``"overview"`` of the lattice, ``"beam"`` for a
        close-up of the beam, ``"racecar"`` for a chase camera in the beam pipe,
        behind the bunch and a little above and to the side of it, looking down the
        pipe ahead, or Ascent camera parameters.
    camera_smoothing : str or float
        For ``projection="3D"`` with ``camera="beam"`` or ``"racecar"``: the width in
        s [m] of a Gaussian average of the reference orbit that the camera follows, so
        that it turns steadily through bends: ``"auto"`` (8 element sizes) or ``0``
        for the local direction of the orbit.
    image_size : tuple of int
        For ``projection="3D"``: width and height of the image in pixels.
    lighting : str, optional
        For ``projection="3D"``: ``None`` draws the translucent elements in their
        colors, ``"camera"`` lights them from the camera, and ``"cinematic"`` lights
        them with a soft sky light, a key light, a fill light and a spot light on the
        beam, in front of a dark background.
    bore : str, optional
        For ``projection="3D"``: ``"round"`` or ``"square"`` bore through the elements,
        or ``None`` for solid elements.
    cutaway : bool
        For ``projection="3D"``: draw the lower half of the elements, cut open through
        the beam axis, opaque. With a bore, drifts are drawn as cut open beam pipes
        around the same bore.
    opacity : float
        For ``projection="3D"``: opacity of a wall of the translucent elements, in
        [0, 1]. ``1`` draws opaque elements, closed around the bore; with a bore,
        drifts are drawn as beam pipes around the same bore. With
        ``camera="racecar"``, they are a tunnel.
    supersample : int
        For ``projection="3D"``: anti-aliasing, the image is rendered this many times
        larger in width and in height and averaged down; ``1`` for none. The rendering
        time grows with its square.
    options : dict, optional
        For ``projection="3D"``: Ascent options, see :py:class:`impactx.insitu.Ascent`.

    Returns
    -------
    Either populates the matplotlib axes in ax or creates a new axes containing the plot.
    With ``projection="3D"``, None on MPI ranks other than 0.

    Raises
    ------
    ValueError
        If the projection is unknown, or for ``projection="3D"``, without a beam and
        without a reference particle, or with both.
    RuntimeError
        For ``projection="3D"``, if ImpactX is built without Ascent, or the envelope
        cannot be tracked through an element.
    """
    if projection == "3D":
        return _plot_survey_3d(
            self,
            ref=ref,
            ax=ax,
            legend=legend,
            legend_ncols=legend_ncols,
            palette=palette,
            beam=beam,
            element=element,
            envelope=envelope,
            magnify=magnify,
            camera=camera,
            camera_smoothing=camera_smoothing,
            image_size=image_size,
            lighting=lighting,
            bore=bore,
            cutaway=cutaway,
            opacity=opacity,
            supersample=supersample,
            options=options,
        )
    if projection != "s":
        raise ValueError(
            f"plot_survey: unknown projection '{projection}', choose 's' or '3D'."
        )

    from math import copysign

    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Rectangle

    from .ElementColors import get_element_color

    charge_qe = 1.0 if ref is None else ref.charge_qe

    ax = ax or plt.subplot(111)

    element_lengths = [element.ds for element in self]

    # NumPy 2.1+ (i.e. Python 3.10+):
    # element_s = np.cumulative_sum(element_lengths, include_initial=True)
    # backport:
    element_s = np.insert(np.cumsum(element_lengths), 0, 0)

    ax.hlines(0, 0, element_s[-1], color="black", linestyle="--")

    # plot config
    skip_names = [
        "Drift",
        "ChrDrift",
        "ExactDrift",
        "Empty",
        "Marker",
        "Source",
    ]

    handles = {}

    for i, element in enumerate(self):
        el_dict = element.to_dict()
        el_type = el_dict["type"]
        el_ds = el_dict["ds"]
        if el_type in skip_names:
            continue

        color = get_element_color(el_type, palette=palette)

        y0 = 0  # default start in y for unspecified elements
        height = 0.5  # default height for unspecified elements

        # note the sub-string matching for el_type
        if el_type == "BeamMonitor":
            y0 = -0.5
            height = 1.0
        if "Quad" in el_type:
            height = copysign(0.8, el_dict["k"] * charge_qe)
        if "Sbend" in el_type:
            if ref is not None:
                height = copysign(0.8, element.rc(ref))
            else:  # guess
                if el_type == "Sbend":
                    el_dict["phi"] = (
                        el_ds / (2 * np.pi * el_dict["rc"]) * 360
                    )  # calculate bending angle (in degrees) and add to dict
                height = copysign(0.8, el_dict["phi"])
        # TODO: sign dependent, read m_p_scale
        # if el_type == "Kicker":
        #    height = copysign(0.8, el_dict["xkick"])

        # plot thin elements on top of thick elements
        zorder = 2
        if element.ds == 0:
            zorder = 3

        patch = Rectangle(
            (element_s[i], y0),
            element_lengths[i],
            height,
            color=color,
            alpha=0.8,
            zorder=zorder,
        )
        ax.add_patch(patch)

        handles[el_type] = patch

    if legend:
        labels = list(handles.keys())
        values = list(handles.values())
        ax.legend(
            handles=values,
            labels=labels,
            bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
            loc="lower left",
            ncols=legend_ncols,
            mode="expand",
            borderaxespad=0.0,
        )

    ax.set_xlabel(r"$s$ [m]")

    ax.set_ylim(-1, 1)
    ax.set_yticks([])

    ax.set_box_aspect(1 / 6)  # some nice aspect ratio

    return ax


def _plot_survey_3d(
    lattice,
    *,
    ref,
    ax,
    legend,
    legend_ncols,
    palette,
    beam,
    element,
    envelope,
    magnify,
    camera,
    camera_smoothing,
    image_size,
    lighting,
    bore,
    cutaway,
    opacity,
    supersample,
    options,
):
    """plot_survey(projection="3D"): Ascent renders the lattice and the beam"""
    import matplotlib.pyplot as plt

    from ..extensions.InSitu import _draw_survey3d, _my_rank, _Renderer
    from ..extensions.InSituSurvey3D import _Survey3D

    vis = _Survey3D(
        lattice,
        ref=ref,
        envelope=envelope,
        magnify=magnify,
        camera=camera,
        image_size=image_size,
        palette=palette,
        lighting=lighting,
        bore=bore,
        cutaway=cutaway,
        opacity=opacity,
        supersample=supersample,
        camera_smoothing=camera_smoothing,
    )
    renderer = _Renderer(options)
    try:
        actions = vis.actions(beam, renderer.image_dir, element=element)
        info = renderer.execute(beam, actions, 0, 0, vis.geometry())
        frame = vis.frame(info, renderer, "survey3d")
    finally:
        renderer.close()

    if _my_rank() != 0:
        return None
    if ax is None:
        width, height = image_size
        dpi = 100
        fig = plt.figure(figsize=(width / dpi, 1.15 * height / dpi), dpi=dpi)
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0 / 1.15])
    _draw_survey3d(
        ax,
        frame["images"][0],
        frame["survey3d"],
        legend=legend,
        legend_ncols=legend_ncols,
    )
    return ax
