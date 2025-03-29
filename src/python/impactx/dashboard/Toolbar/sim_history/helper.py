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

