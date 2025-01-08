from .. import TrameFunctions, generalFunctions, setup_server, vuetify

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------


@state.change("csr_bins")
def on_csr_bins_change(csr_bins, **kwargs):
    state.csr_bins_error_message = generalFunctions.validate_against(
        csr_bins, "int", ["positive"]
    )
    generalFunctions.update_simulation_validation_status()


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class csrConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for CSR.
        """

        with vuetify.VCard(v_show="csr", style="width: 170px;"):
            TrameFunctions.input_section_header("CSR")
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        TrameFunctions.select(
                            label="Particle Shape",
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        TrameFunctions.text_field(
                            label="CSR Bins",
                        )
