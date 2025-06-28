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
    def select_visualization() -> vuetify.VTabs:
        """
        Provides the user a tab group to select the type of visualization to view.
        """
        return vuetify.VTabs(
            v_model=("active_visualization",),
            items=("visualization_options",),
            color="primary",
            hide_slider=False,
            disabled=("!sims.length",),  # disabled if no sims are in the history
        )
