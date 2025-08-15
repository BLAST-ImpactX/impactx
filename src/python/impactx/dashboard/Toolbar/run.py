"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import html

from .. import ctrl, state, vuetify
from ..Input.components import CardComponents
from ..Run.executor import run_execute_impactx_sim


class RunToolbar:
    """
    Contains toolbar components for the 'Run' page.
    """

    @ctrl.trigger("begin_sim")
    def run():
        """
        Called when the 'Run Simulation' button is clicked.
        """
        run_execute_impactx_sim()

    @ctrl.trigger("cancel_sim")
    def cancel_sim():
        """
        Called when the 'Run Simulation' button is clicked while
        a simulation is on-going.
        """
        state.sim_is_cancelled = True

    @staticmethod
    def run_simulation():
        """
        Displays the 'Run Simulation' components
        """
        (RunToolbar.run_simulation_progress_details(),)
        (RunToolbar.run_simulation_progress_bar(),)
        (RunToolbar.run_simulation_button(),)

    @staticmethod
    def run_simulation_button() -> vuetify.VBtn:
        """
        Component to run the simulation.

        On click, it either starts the simulation
        or cancels it if it is already running.

        Disabled when the simulation is generating plots
        or if the simulation is not running.

        Color changes based on the simulation status.
        """

        CardComponents.card_button(
            ["mdi-play-circle", "mdi-close-circle"],
            color=("sim_is_running ? 'error' : sim_status_color",),
            click="sim_is_running ? trigger('cancel_sim') : trigger('begin_sim')",
            description=["Run Simulation", "Cancel Simulation"],
            dynamic_condition="sim_is_running",
            disabled=("disableRunSimulationButton || sim_is_generating_plots", True),
            id="run_simulation_button",
        )

    @staticmethod
    def run_simulation_progress_bar() -> vuetify.VProgressLinear:
        """
        Displays a progress bar indicating the current progress of the simulation.

        Below the progress bar, it shows the current simulation status
        (e.g., "Running", "Completed", "Cancelled") and the percentage of completion.
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
        Displays the current step and elapsed time of the simulation.
        This component is updated dynamically to reflect the current state of the simulation.
        """

        return html.Div(
            "Step {{ sim_current_step }} • {{ sim_elapsed_time }}",
            style="margin-right: 8px;",
        )
