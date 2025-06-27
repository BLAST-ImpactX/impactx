"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from .. import html, setup_server, vuetify

server, state, ctrl = setup_server()

class AnalyzeToolbar:
    """
    Contains toolbar components for the 'Analyze' page.
    """

    @staticmethod
    def select_visualization() -> vuetify.VSelect:
        """
        Provides the user a dropdown to select the type of visualization to view.
        """
        with html.Div(style="width: 15vw"):
            vuetify.VSelect(
                v_model=("active_visualization",),
                items=("visualization_options",),
                label="Select Visualization",
                disabled=("!sims.length",),  # disabled if no sims are in the history
                density="comfortable",
                variant="solo-filled",
                hide_details=True,
                color="primary",
            )
