"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import html, setup_server, vuetify
from .components import SimulationHistoryComponents
from ...Input.components.navigation import NavigationComponents
server, state, ctrl = setup_server()


@staticmethod
def sim_details_tabs():
    dialog_name = "sim_details_tabs"
    with NavigationComponents.create_dialog_tabs(dialog_name, 1, ["Inputs"]):
        with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
            with vuetify.VTabsWindowItem():
                with vuetify.VCardText():
                    with html.Div(classes="code-editor-style"):
                        html.Div("{{ selected_sim?.inputs }}")
                    
class SimulationHistoryDialogs:

    @staticmethod
    def rename_dialog():
        """
        Contains the UI and functionality for the
        simulation history 'Rename' action button.
        """
        with vuetify.VDialog(v_model=("sim_rename_dialog", False), max_width=400, persistent=True):
            with vuetify.VCard():
                vuetify.VCardTitle("Rename Simulation")
                with vuetify.VCardText():
                    with vuetify.VRow():
                        with vuetify.VCol():
                            SimulationHistoryComponents.text_field(
                                label="Current Name",
                                v_model_name="rename_old_name",
                                readonly=True,
                                disabled=True,
                            )
                    with vuetify.VRow():
                        with vuetify.VCol():
                            SimulationHistoryComponents.text_field(
                                label="New Name",
                                v_model_name="rename_new_name",
                                clearable=True,
                            )
                with vuetify.VCardActions():
                    vuetify.VSpacer()
                    vuetify.VBtn(
                        "Cancel",
                        variant="text",
                        click=ctrl.close_rename_dialog
                    )
                    vuetify.VBtn(
                        "Confirm",
                        color="primary",
                        variant="elevated",
                        click=ctrl.confirm_rename,
                    )

    @staticmethod
    def sim_details_dialog():
        """
        Contains the UI and functionality for the
        simulation history 'View Details' action button.
        """
        with vuetify.VDialog(v_model=("sim_details_dialog", False), max_width="700px"):
            with vuetify.VCard(elevation=10, classes="rounded-lg"):
                with vuetify.VToolbar(color="primary", classes="px-4"):
                    vuetify.VIcon("mdi-clipboard-text-clock")
                    vuetify.VToolbarTitle("{{ selected_sim?.name }}")
                    vuetify.VSpacer()
                    vuetify.VBtn(icon="mdi-close", click="sim_details_dialog = false")
                with vuetify.VCardText():
                    with html.Div(classes="ga-4 d-flex flex-wrap mb-2"):
                        with SimulationHistoryComponents.sim_details_card(
                            title="STATUS"
                        ):
                            with html.Div():
                                SimulationHistoryComponents.status_chip("selected_sim?")
                        with SimulationHistoryComponents.sim_details_card(
                            title="CREATED",
                            prepend_icon="mdi-calendar"
                        ):
                            with html.Div():
                                html.Span("{{ new Date(selected_sim?.created_at_time).toLocaleString() }}", classes="font-weight-medium")
                        with SimulationHistoryComponents.sim_details_card(title="DURATION", prepend_icon="mdi-clock-outline"):
                            with html.Div():
                                html.Span(
                                    "{{ selected_sim?.time_elapsed || '—' }}",
                                    classes="font-weight-medium"
                            )
                    with vuetify.VCard(elevation=2):
                        sim_details_tabs()
