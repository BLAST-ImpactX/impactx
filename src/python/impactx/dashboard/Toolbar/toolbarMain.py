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
from .importParser import DashboardParser

server, state, ctrl = setup_server()

state.show_dashboard_alert = True
state.import_file = False
state.import_file_details = None
state.import_file_error = False
state.importing_file = False


# -----------------------------------------------------------------------------
# Triggers/Controllers
# -----------------------------------------------------------------------------
def reset_importing_states():
    state.import_file_error = None
    state.import_file_details = None
    state.import_file = None
    state.importing_file = False


@ctrl.add("reset_all")
def reset_all():
    reset_importing_states()
    generalFunctions.reset_inputs("all")


@ctrl.trigger("export")
def on_export_click():
    return input_file()


@state.change("import_file")
def on_import_file_change(import_file, **kwargs):
    if import_file:
        try:
            state.importing_file = True
            DashboardParser.file_details(import_file)
            DashboardParser.populate_impactx_simulation_file_to_ui(import_file)
        except Exception:
            state.import_file_error = True
            state.import_file_error_message = "Unable to parse"
        finally:
            state.importing_file = False


# -----------------------------------------------------------------------------
# Common toolbar elements
# -----------------------------------------------------------------------------


class ToolbarElements:
    """
    Helper functions to create
    Vuetify UI elements for toolbar.
    """

    @staticmethod
    def export_button():
        with vuetify.VBtn(
            click="utils.download('impactx_simulation.py', trigger('export'), 'text/plain')",
            outlined=True,
            small=True,
            disabled=("disableRunSimulationButton", True),
            classes="mx-2",
        ):
            vuetify.VIcon("mdi-download", left=True, small=True)
            html.Span("Export")

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
    def import_button():
        vuetify.VFileInput(
            v_model=("import_file",),
            accept=".py",
            __properties=["accept"],
            style="display: none;",
            ref="fileInput",
        )
        with html.Div(
            style="position: relative",
        ):
            with vuetify.VBtn(
                click="$refs.fileInput.$refs.input.click()",
                outlined=True,
                small=True,
                disabled=("(import_file_details)",),
                color=("import_file_error ? 'error' : ''",),
            ):
                vuetify.VIcon(
                    "mdi-upload",
                    left=True,
                    small=True,
                )
                html.Span("Import")
            with html.Div(
                style="position: absolute; font-size: 10px; width: 100%; padding-top: 2px; display: flex; justify-content: center; white-space: nowrap;"
            ):
                html.Span(
                    "{{ import_file_error ? import_file_error_message : import_file_details }}",
                    style="text-overflow: ellipsis; overflow: hidden;",
                    classes=(
                        "import_file_error ? 'error--text' : 'grey--text text--darken-1'",
                    ),
                )
                vuetify.VIcon(
                    "mdi-close",
                    x_small=True,
                    style="cursor: pointer;",
                    click=ctrl.reset_all,
                    v_if="import_file_details || import_file_error",
                    color=("import_file_error ? 'error' : 'grey darken-1'",),
                )

    @staticmethod
    def reset_inputs_button():
        with vuetify.VBtn(
            click=ctrl.reset_all,
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
            ToolbarElements.import_button()
            ToolbarElements.export_button()
            ToolbarElements.reset_inputs_button()

        elif toolbar_name == "run":
            (ToolbarElements.dashboard_info(),)
            (vuetify.VSpacer(),)
            (ToolbarElements.run_simulation_button(),)
        elif toolbar_name == "analyze":
            (ToolbarElements.dashboard_info(),)
            vuetify.VSpacer()
            ToolbarElements.plot_options()
