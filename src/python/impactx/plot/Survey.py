"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def plot_survey(
    self, ref=None, ax=None, legend=True, legend_ncols=5, palette="cern-lhc"
):
    """Plot over s of all elements in the KnownElementsList.

    A positive element strength denotes horizontal focusing (e.g. for quadrupoles) and bending to the right (for dipoles).  In general, this depends on both the sign of the field and the sign of the charge.

    Parameters
    ----------
    self : ImpactXParticleContainer_*
        The KnownElementsList class in ImpactX
    ref : RefPart
        A reference particle, checked for the charge sign to plot focusing/defocusing strength directions properly.
    ax : matplotlib axes
        A plotting area in matplotlib (called axes there).
    legend: bool
        Plot a legend if true.
    legend_ncols: int
        Number of columns for lattice element types in the legend.
    palette: string
        Color palette.

    Returns
    -------
    Either populates the matplotlib axes in ax or creates a new axes containing the plot.
    """
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
            if ref is None:
                height = copysign(0.8, element.rc(ref))
            else:  # guess
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

    ax.set_aspect(1 / 1.618)  # golden ratio

    return ax
