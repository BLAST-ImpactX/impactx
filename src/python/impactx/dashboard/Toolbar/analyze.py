"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from .. import setup_server, vuetify

server, state, ctrl = setup_server()


class AnalyzeToolbar:
    """
    Contains toolbar components for the 'Analyze' page.
    """

    @staticmethod
    def plot_options() -> vuetify.VSelect:
        """
        Displays a dropdown menu to select the available plots
        generated from the simulation data.
        """

        return vuetify.VSelect(
            v_model=("active_plot", "1D plots over s"),
            items=("plot_options",),
            label="Select plot to view",
            hide_details=True,
            density="compact",
            variant="underlined",
            style="max-width: 250px",
            disabled=("disableRunSimulationButton", True),
        )
