"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import html, vuetify

from ..Input.generalFunctions import generalFunctions
from ..trame_setup import setup_server
from .exportTemplate import input_file

server, state, ctrl = setup_server()

state.show_dashboard_alert = True

# -----------------------------------------------------------------------------
# Triggers
# -----------------------------------------------------------------------------


@ctrl.trigger("export")
def on_export_click():
    return input_file()


# -----------------------------------------------------------------------------
# Common toolbar elements
# -----------------------------------------------------------------------------


class ToolbarElements:
    """
    Helper functions to create
    Vuetify UI elements for toolbar.
    """

    @staticmethod
    def export_input_data():
        vuetify.VIcon(
            "mdi-download",
            style="color: #00313C;",
            click="utils.download('impactx_simulation.py', trigger('export'), 'text/plain')",
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def plot_options():
        vuetify.VSelect(
            v_model=("active_plot", "1D plots over s"),
            items=("plot_options",),
            label="Select plot to view",
            hide_details=True,
            dense=True,
            style="max-width: 250px",
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def run_simulation_button():
        vuetify.VBtn(
            "Run Simulation",
            style="background-color: #00313C; color: white; margin: 0 20px;",
            click=ctrl.run_simulation,
            disabled=("disableRunSimulationButton", True),
        )

    @staticmethod
    def reset_inputs_button():
        with vuetify.VBtn(
            click=lambda: generalFunctions.reset_inputs("all"),
            outlined=True,
            small=True,
        ):
            vuetify.VIcon("mdi-refresh", left=True)
            html.Span("Reset")

    @staticmethod
    def dashboard_info():
        """
        Creates an alert box with dashboard information.
        """

        vuetify.VAlert(
            "ImpactX Dashboard is provided as a preview and continues to be developed. "
            "Thus, it may not yet include all the features available in ImpactX.",
            type="info",
            dense=True,
            dismissible=True,
            v_model=("show_dashboard_alert", True),
            classes="mt-4",
        )


class Toolbars:
    """
    Builds toolbar for dashboard.
    """

    @staticmethod
    def dashboard_toolbar(toolbar_name: str) -> None:
        toolbar_name = toolbar_name.lower()
        if toolbar_name == "input":
            (ToolbarElements.dashboard_info(),)
            vuetify.VSpacer()
            ToolbarElements.reset_inputs_button()
            ToolbarElements.export_input_data()
        elif toolbar_name == "run":
            (ToolbarElements.dashboard_info(),)
            (vuetify.VSpacer(),)
            (ToolbarElements.run_simulation_button(),)
        elif toolbar_name == "analyze":
            (ToolbarElements.dashboard_info(),)
            vuetify.VSpacer()
            ToolbarElements.plot_options()
