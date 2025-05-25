"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import html

from .. import setup_server, vuetify
from ..Analyze.plotsMain import available_plot_options, load_dataTable_data, update_plot
from ..Input.components.card import CardComponents
from ..Run.executor import run_execute_impactx_sim
from .dashboard_parser.parser import DashboardParser
from .sim_history.ui import SimulationHistory

server, state, ctrl = setup_server()

state.show_dashboard_alert = True
state.import_file = False
state.import_file_details = None
state.import_file_error = False
state.importing_file = False

state.expand_all_sections = False



from .input_toolbar import InputToolbar

class RunToolbar:
    """
    Contains toolbar elements for the Run page.
    """

    @ctrl.trigger("begin_sim")
    def run():
        state.plot_options = available_plot_options(simulationClicked=True)
        run_execute_impactx_sim()
        update_plot()
        load_dataTable_data()

    @ctrl.trigger("cancel_sim")
    def cancel_sim():
        state.sim_is_cancelled = True

    @staticmethod
    def run_simulation():
        (RunToolbar.run_simulation_progress_details(),)
        (RunToolbar.run_simulation_progress_bar(),)
        (RunToolbar.run_simulation_button(),)

    @staticmethod
    def run_simulation_button() -> vuetify.VBtn:
        """
        Creates a button to run an ImpactX simulation
        with the current user-provided inputs.
        """
        CardComponents.card_button(
            ["mdi-play-circle", "mdi-close-circle"],
            color=("sim_is_running ? 'error' : sim_status_color",),
            click="sim_is_running ? trigger('cancel_sim') : trigger('begin_sim')",
            description="Run Simulation",
            dynamic_condition="sim_is_running",
            disabled=("disableRunSimulationButton || sim_is_generating_plots", True),
        )

    @staticmethod
    def run_simulation_progress_bar() -> vuetify.VBtn:
        """
        Displays and updates a progress bar to the dashboard user
        while running a simulation.
        """
        with html.Div(style="position: relative; margin: 0 8px;"):
            vuetify.VProgressLinear(
                height=5,
                striped=True,
                style="width: 7vw",
                color=("sim_status_color",),
                v_model=("sim_progress",),
            )
            html.Div(
                "{{ sim_progress_status }}",
                style="position: absolute; top: 100%; left: 50%; transform: translateX(-50%); font-size: 12px; white-space: nowrap; color: grey; margin-top: 4px;",
            )

    @staticmethod
    def run_simulation_progress_details() -> html.Div:
        """
        Provides dashboard users with simulation progress details,
        such as the current step and the time elapsed in the simulation.
        """

        return html.Div(
            "Step {{ sim_current_step }} • {{ sim_elapsed_time }}",
            style="margin-right: 8px;",
        )


class AnalyzeToolbar:
    """
    Contains toolbar elements for the Analyze page.
    """

    @staticmethod
    def plot_options() -> vuetify.VSelect:
        """
        Creates a dropdown menu for selecting a plot
        to visualize simulation results.
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


class GeneralToolbar:
    """
    General tolbar elements.
    """

    @staticmethod
    def dashboard_toolbar(toolbar_name: str) -> None:
        """
        Builds and displays the appropriate toolbar
        based on the selected dashboard section.

        :param toolbar_name: The name of the dashboard section
        for which the toolbar is needed.
        """

        toolbar_name = toolbar_name.lower()
        if toolbar_name == "input":
            (GeneralToolbar.dashboard_info(),)
            vuetify.VSpacer()
            InputToolbar.import_button()
            InputToolbar.export_button()
            InputToolbar.reset_inputs_button()
            vuetify.VDivider(vertical=True, classes="mr-2")
            GeneralToolbar.simulation_history_button()
            vuetify.VDivider(vertical=True, classes="mr-2")
            InputToolbar.collapse_all_sections_button()
        elif toolbar_name == "run":
            (GeneralToolbar.dashboard_info(),)
            (vuetify.VSpacer(),)
            (RunToolbar.run_simulation(),)
            vuetify.VDivider(vertical=True, classes="mx-2")
            (GeneralToolbar.simulation_history_button())
        elif toolbar_name == "analyze":
            (GeneralToolbar.dashboard_info(),)
            vuetify.VSpacer()
            AnalyzeToolbar.plot_options()

    @staticmethod
    def dashboard_info() -> vuetify.VAlert:
        """
        Creates an informational alert box for the dashboard to
        notify users that the ImpactX dashboard is still in development.

        :return: A Vuetify alert component displaying the dashboard notice.
        """

        return vuetify.VAlert(
            "ImpactX Dashboard is provided as a preview and continues to be developed. "
            "Thus, it may not yet include all the features available in ImpactX.",
            type="info",
            density="compact",
            dismissible=True,
            v_model=("show_dashboard_alert", True),
            classes="text-body-2 hidden-md-and-down",
            style="width: 50vw; overflow: hidden; margin: auto;",
        )

    @staticmethod
    def simulation_history_button() -> vuetify.VBtn:
        """
        Displays a button to the user which holds the
        components to the simulation history.
        """
        SimulationHistory.simulation_history()
        SimulationHistory.init_sim_history_dialogs()

        return vuetify.VBtn(
            "History",
            color="primary",
            classes="mr-2",
            click="simulation_history_dialog = true",
            prepend_icon="mdi-clipboard-text-clock",
            size="small",
            variant="elevated",
            disabled=("!sims.length",),
        )
