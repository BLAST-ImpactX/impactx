"""

This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

__all__: list[str] = ["plot_survey"]

def plot_survey(
    self, ref=None, ax=None, legend=True, legend_ncols=5, palette="cern-lhc"
):
    """
    Plot over s of all elements in the KnownElementsList.

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
