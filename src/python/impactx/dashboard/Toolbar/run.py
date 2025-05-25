"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import html

from .. import setup_server, vuetify
from ..Analyze.plotsMain import available_plot_options, load_dataTable_data, update_plot
from ..Input.components.card import CardComponents
from ..Run.executor import run_execute_impactx_sim

server, state, ctrl = setup_server()

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