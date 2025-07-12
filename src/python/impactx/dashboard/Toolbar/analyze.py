"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from .. import state, vuetify

state.selected_sim_to_analyze = None


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

    @staticmethod
    def simulation_selection_indicator() -> vuetify.VChip:
        """
        Displays the selected simulation for analysis.

        By default, it shows the most recently run simulation if one is available.
        """

        return vuetify.VChip(
            "{{ sim_is_running ? sim_progress_status : (selected_sim_to_analyze?.name || 'No simulation') }}",
            color=("sim_is_running ? 'info' : 'green-darken-1'",),
            prepend_icon="mdi-check-circle-outline",
        )
