from .. import TrameFunctions, setup_server, vuetify

server, state, ctrl = setup_server()


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
                            input=(ctrl.input_change, "['csr_bins']"),
                        )
