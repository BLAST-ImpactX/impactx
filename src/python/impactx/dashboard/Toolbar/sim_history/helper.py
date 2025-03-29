"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import html, setup_server, vuetify

server, state, ctrl = setup_server()

class SimulationHistoryDialogs:

    @staticmethod
    def rename_dialog():
        with vuetify.VDialog(v_model=("sim_rename_dialog", False), max_width=400, persistent=True):
            with vuetify.VCard():
                vuetify.VCardTitle("Rename Simulation")
                with vuetify.VCardText():
                    with vuetify.VRow():
                        with vuetify.VCol():
                            vuetify.VTextField(
                                label="Current Name",
                                v_model=("rename_old_name", ""),
                                readonly=True,
                                hide_details=True,
                                variant="outlined",
                                density="comfortable",
                                disabled=True,
                            )
                    with vuetify.VRow():
                        with vuetify.VCol():
                            vuetify.VTextField(
                                label="New Name",
                                v_model=("rename_new_name", ""),
                                hide_details=True,
                                variant="outlined",
                                density="comfortable",
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

    def sim_details_dialog():
        with vuetify.VDialog(v_model=("sim_details_dialog", False), style="width: 500px"):
            with vuetify.VCard(elevation=6):
                with vuetify.VToolbar(color="primary", classes="px-4"):
                    vuetify.VToolbarTitle("Simulation Details")
                    vuetify.VSpacer()
                    vuetify.VBtn(
                        icon="mdi-close",
                        click="sim_details_dialog = false",
                    )
                with vuetify.VCardText():
                    html.P("Name: {{ selected_sim?.name }}")
                    html.P("Status: {{ selected_sim?.status }}")
                    html.P("Created At: {{ new Date(selected_sim?.created_at_time).toLocaleString() }}")
                    html.P("Duration: {{ selected_sim?.time_elapsed }}")