"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import plotly

from .. import ctrl, state, vuetify
from ..Input.components.navigation import NavigationComponents
from .over_s.selector import over_s_selector

state.visualization_options = ["Plot Over S", "Phase Space Plots"]
state.active_visualization = "Plot Over S"
state.phase_space_png = None


class AnalyzeSimulation:
    """
    Prepares contents for the 'Analyze' page.
    """

    @staticmethod
    def plot_over_s():
        """
        Displays the content for the 'Plot Over S' selection.
        """

        dialog_name = "plot_over_s_tab_dialog"

        with vuetify.VContainer(fluid=True):
            with vuetify.VRow():
                with vuetify.VCol(cols=9, classes="d-flex flex-column"):
                    with NavigationComponents.create_dialog_tabs(
                        dialog_name, 2, ["Plot", "Data"]
                    ):
                        with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
                            with vuetify.VTabsWindowItem():  # tab1
                                with vuetify.VContainer(
                                    style="height: calc(84vh - 8px); width: 100%;",
                                ):
                                    plotly_figure = plotly.Figure(
                                        display_mode_bar="true",
                                    )
                                    ctrl.plotly_figure_update = plotly_figure.update
                            with vuetify.VTabsWindowItem():  # tab2
                                with vuetify.VContainer(
                                    style="height: calc(84vh - 8px); width: 100%;",
                                ):
                                    vuetify.VDataTable(
                                        headers=("over_s_table_headers",),
                                        items=("over_s_table_data", []),
                                        density="compact",
                                    )
                with vuetify.VCol(cols=3, classes="fill-height"):
                    over_s_selector.selector()

    @staticmethod
    def phase_space():
        """
        Displays the phase space plots.
        """

        with vuetify.VContainer(fluid=True):
            with vuetify.VContainer():
                with vuetify.VCard(style="height: 50vh; width: 150vh;"):
                    with vuetify.VTabs(v_model=("active_tab", 0)):
                        vuetify.VTab("Plot")
                    vuetify.VDivider()
                    with vuetify.VTabsWindow(v_model="active_tab"):
                        with vuetify.VTabsWindowItem():
                            vuetify.VImg(
                                v_if=("phase_space_png",), src=("phase_space_png",)
                            )
