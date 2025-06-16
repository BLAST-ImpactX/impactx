"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""


from trame.widgets import plotly

from . import over_s
from .. import setup_server, vuetify
from ..Input.components.navigation import NavigationComponents

server, state, ctrl = setup_server()

PLOTS = {
    "Plot Over S": over_s.update,
    "Phase Space Plots": None,
}


def available_plot_options(simulationClicked):
    """
    Displays plot_options for users based on status of simulation.
    :param simulationClicked (bool): status of simulation status
    :return: list of plot_options for users
    """

    if simulationClicked:
        return list(PLOTS.keys())
    else:
        return ["Run Simulation To See Options"]


state.plot_options = available_plot_options(simulationClicked=False)
state.active_plot = None
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
                    with NavigationComponents.create_dialog_tabs(dialog_name, 2, ["Plot", "Data"]):
                        with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
                            with vuetify.VTabsWindowItem(): # tab1
                                with vuetify.VContainer(style="height: 80vh; width: 100%;",):
                                    plotly_figure = plotly.Figure(
                                        display_mode_bar="true",
                                    )
                                    ctrl.plotly_figure_update = plotly_figure.update
                            with vuetify.VTabsWindowItem(): # tab2
                                with vuetify.VContainer(style="height: 80vh; width: 100%;",):
                                    vuetify.VDataTable(
                                        headers=("over_s_table_headers",),
                                        items=("over_s_table_data", []),
                                        density="compact",
                                    )
                with vuetify.VCol(cols=3):
                    with vuetify.VCard(classes="pa-4 d-flex flex-column"):
                        vuetify.VSelect(
                            v_model=("selected_headers",),
                            items=("selectable_headers",),
                            label="Select data to view",
                            multiple=True,
                            density="comfortable",
                            variant="outlined",
                        )

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
