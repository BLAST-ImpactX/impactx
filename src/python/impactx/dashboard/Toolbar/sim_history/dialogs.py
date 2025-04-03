"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import html, setup_server, vuetify
from ...Input.components.navigation import NavigationComponents
from .components import SimulationHistoryComponents

server, state, ctrl = setup_server()


@staticmethod
def view_details_tabs():
    dialog_name = "view_details_tabs"
    with NavigationComponents.create_dialog_tabs(
        dialog_name, 2, ["Log", "Inputs File"]
    ):
        with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
            with vuetify.VTabsWindowItem():
                with vuetify.VCardText():
                    with html.Div(classes="code-editor-style"):
                        html.Div("{{ selected_sim?.log }}")
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
        with SimulationHistoryComponents.dialog(
            title="{{ selected_sim?.name }} - Rename",
            prepend_icon="mdi-pencil",
            dialog_var="sim_rename_dialog",
            width="33.33vw",
        ):
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
                    "Confirm",
                    color="primary",
                    variant="elevated",
                    click=ctrl.confirm_rename,
                )

    @staticmethod
    def view_details_dialog():
        """
        Contains the UI and functionality for the
        simulation history 'View Details' action button.
        """

        with SimulationHistoryComponents.dialog(
            title="{{ selected_sim?.name }} - Details",
            prepend_icon="mdi-clipboard-text-clock",
            dialog_var="view_details_dialog",
            width="70vw",
        ):
            with vuetify.VCardText():
                with html.Div(classes="ga-4 d-flex flex-wrap mb-2"):
                    with SimulationHistoryComponents.view_details_card(
                        title="STATUS", prepend_icon="mdi-label-outline"
                    ):
                        with html.Div():
                            SimulationHistoryComponents.status_chip("selected_sim?")
                    with SimulationHistoryComponents.view_details_card(
                        title="CREATED", prepend_icon="mdi-calendar"
                    ):
                        with html.Div():
                            html.Span(
                                "{{ window.formatDate(selected_sim.created_at_time) }} at {{ window.formatTime(selected_sim.created_at_time) }}",
                                classes="font-weight-medium",
                            )
                    with SimulationHistoryComponents.view_details_card(
                        title="DURATION", prepend_icon="mdi-clock-outline"
                    ):
                        with html.Div():
                            html.Span(
                                "{{ selected_sim?.time_elapsed || '—' }}",
                                classes="font-weight-medium",
                            )
                with vuetify.VCard(elevation=2):
                    view_details_tabs()

    @staticmethod
    def download_options_dialog():
        with SimulationHistoryComponents.dialog(
            title="{{ sim_to_download?.name }} - Downloading Options",
            prepend_icon="mdi-download",
            dialog_var="sim_download_dialog",
            width="33.33vw",
        ):
            with vuetify.VCardText():
                with vuetify.VList():
                    vuetify.VListItem(
                        title="Download Inputs",
                        prepend_icon="mdi-file-code",
                        click="utils.download(`${sim_to_download.name}.py`, trigger('download_sim', [sim_to_download]), 'text/plain'); sim_download_dialog = false",
                    )

    @staticmethod
    def load_sim_dialog():
        with SimulationHistoryComponents.dialog(
            title="{{ selected_sim_to_load?.name }} - Loading Options",
            prepend_icon="mdi-upload",
            dialog_var="load_sim_dialog",
            width="33.33vw",
        ):
            with vuetify.VCardText():
                with vuetify.VList():
                    vuetify.VListItem(
                        title="Load Inputs",
                        prepend_icon="mdi-file-code",
                        click=(ctrl.load_selected_sim),
                    )
