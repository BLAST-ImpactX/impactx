"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import os

from .. import ctrl, html, state, vuetify
from ..Input.components import CardComponents
from .analyze import AnalyzeToolbar
from .input import InputToolbar
from .run import RunToolbar
from .sim_history.ui import SimulationHistory

state.show_dashboard_alert = True
state.about_dialog = False


@ctrl.trigger("force_quit")
def _force_quit() -> None:
    os._exit(0)


def _about_button(
    icon_name: str, text: str, link: str, color: str = "primary"
) -> vuetify.VBtn:
    """
    Returns a button for the about section.
    """
    return vuetify.VBtn(
        prepend_icon=icon_name,
        text=text,
        color=color,
        classes="justify-start text-none",
        variant="outlined",
        click=f"window.open('{link}', '_blank')",
    )


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
            GeneralToolbar.about_button()
            GeneralToolbar.force_quit_button()
        elif toolbar_name == "run":
            (GeneralToolbar.dashboard_info(),)
            (vuetify.VSpacer(),)
            (RunToolbar.run_simulation(),)
            vuetify.VDivider(vertical=True, classes="mx-2")
            (GeneralToolbar.simulation_history_button())
            vuetify.VDivider(vertical=True, classes="mx-2")
            GeneralToolbar.about_button()
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
            GeneralToolbar.about_button()
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

    @staticmethod
    def about_button() -> vuetify.VBtn:
        """
        Displays a button to open the about dialog.
        """
        GeneralToolbar.about_dialog()

        return CardComponents.card_button(
            icon_name="mdi-information-outline",
            click="about_dialog = true",
            description="About",
            color="primary",
            classes="mr-2",
        )

    @staticmethod
    def about_dialog() -> None:
        """
        Creates the about dialog with information and links to ImpactX resources.
        """
        HEADER_1 = "About ImpactX Dashboard"
        MESSAGE_1 = "This dashboard is a web-based interface for monitoring and analyzing particle accelerator simulations based on Trame."
        MESSAGE_2 = "This dashboard provides visualization and analysis tools for ImpactX simulations, however it is a subset of the complete ImpactX ecosystem."
        HEADER_2 = "Get the Full ImpactX Experience"
        MESSAGE_3 = "For complete simulation capabilities and local installations:"
        HEADER_3 = "Documentation"
        BUG_MESSAGE = "Found a bug or have a feature request? "

        IMPACTX_DOCUMENTATION_URL = "https://impactx.readthedocs.io/"
        IMPACTX_EXAMPLES_URL = (
            "https://impactx.readthedocs.io/en/latest/usage/examples.html"
        )
        GITHUB_URL = "https://github.com/BLAST-ImpactX/impactx"
        GITHUB_ISSUES_URL = "https://github.com/BLAST-ImpactX/impactx/issues/new"
        DISCUSSIONS_URL = "https://github.com/orgs/BLAST-ImpactX/discussions"

        with vuetify.VDialog(v_model=("about_dialog", False), max_width="500px"):
            with vuetify.VCard(elevation=10, classes="rounded-lg"):
                with vuetify.VToolbar(color="primary", classes="px-4"):
                    html.Div(
                        HEADER_1,
                        style="font-size: 1.125rem;",
                    )
                    vuetify.VSpacer()
                    vuetify.VBtn(icon="mdi-close", click="about_dialog = false")

                with vuetify.VCardText(classes="pa-6"):
                    with html.Div(classes="mb-4"):
                        html.P(MESSAGE_1, classes="mb-3")
                        html.P(MESSAGE_2, classes="mb-4")
                    vuetify.VDivider(classes="my-4")
                    with html.Div(classes="mb-4"):
                        html.H4(HEADER_2, classes="text-h6 mb-3")
                        html.P(MESSAGE_3, classes="mb-3 text-body-2")

                    with html.Div(classes="d-flex flex-column ga-3 mb-4"):
                        _about_button(
                            icon_name="mdi-book-open-variant",
                            text=HEADER_3,
                            link=IMPACTX_DOCUMENTATION_URL,
                        )
                        _about_button(
                            icon_name="mdi-github",
                            text="Source Code",
                            color="secondary",
                            link=GITHUB_URL,
                        )
                        _about_button(
                            icon_name="mdi-play-circle-outline",
                            text="Examples",
                            link=IMPACTX_EXAMPLES_URL,
                            color="info",
                        )
                        _about_button(
                            icon_name="mdi-comment-question-outline",
                            text="Questions & Answers",
                            link=DISCUSSIONS_URL,
                            color="settings",
                        )

                    with vuetify.VAlert(
                        type="info", variant="tonal", density="compact"
                    ):
                        with vuetify.Template(v_slot_prepend=True):
                            vuetify.VIcon("mdi-bug-outline")
                        html.Span(BUG_MESSAGE)
                        html.A(
                            "Report it on GitHub",
                            href=GITHUB_ISSUES_URL,
                            target="_blank",
                            classes="text-decoration-underline",
                            style="color: #1976d2; cursor: pointer;",
                        )
