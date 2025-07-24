"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import os

from .. import ctrl, state, vuetify
from ..Input.components import CardComponents
from .sim_history.ui import SimulationHistory

state.show_dashboard_alert = True

from .analyze import AnalyzeToolbar
from .input import InputToolbar
from .run import RunToolbar


@ctrl.trigger("force_quit")
def _force_quit() -> None:
    os._exit(0)


class GeneralToolbar:
    """
    Contains toolbar components displayed on all pages.
    """

    @staticmethod
    def dashboard_toolbar(toolbar_name: str) -> None:
        """
        Displays the toolbar components based on the provided toolbar name.
        The toolbar name should be one of the following:
        - "input": Displays components related to input configuration.
        - "run": Displays components related to running simulations.
        - "analyze": Displays components related to analyzing simulation results.

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
            GeneralToolbar.force_quit_button()
        elif toolbar_name == "run":
            (GeneralToolbar.dashboard_info(),)
            (vuetify.VSpacer(),)
            (RunToolbar.run_simulation(),)
            vuetify.VDivider(vertical=True, classes="mx-2")
            (GeneralToolbar.simulation_history_button())
            vuetify.VDivider(vertical=True, classes="mx-2")
            (GeneralToolbar.force_quit_button())
        elif toolbar_name == "analyze":
            (GeneralToolbar.dashboard_info(),)
            vuetify.VSpacer()
            AnalyzeToolbar.select_visualization()
            vuetify.VDivider(vertical=True, classes="mx-2")
            AnalyzeToolbar.simulation_selection_indicator()
            vuetify.VDivider(vertical=True, classes="mx-2")
            GeneralToolbar.simulation_history_button()
            vuetify.VDivider(vertical=True, classes="mx-2")
            GeneralToolbar.force_quit_button()

    @staticmethod
    def dashboard_info() -> vuetify.VAlert:
        """
        Displays an informational alert box for the dashboard to
        notify users that the ImpactX dashboard is still in development.
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
        Displays a button to open the simulation history dialog.

        This button is disabled when there are no simulations available
        (ie, when `sims.length` is 0).
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

    @staticmethod
    def force_quit_button():
        """
        Displays a button to force quit the dashboard.
        """
        return CardComponents.card_button(
            icon_name="mdi-power",
            click="trigger('force_quit')",
            description="Force Quit",
            color="error",
        )
