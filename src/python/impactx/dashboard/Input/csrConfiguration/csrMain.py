from .. import CardComponents, InputComponents, setup_server, vuetify

server, state, ctrl = setup_server()


class csrConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for CSR.
        """

        with vuetify.VCard(v_show="csr", style="width: 170px;"):
            CardComponents.input_header("CSR")
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        InputComponents.select(
                            label="Particle Shape",
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(classes="py-0"):
                        InputComponents.text_field(
                            label="CSR Bins",
                            input=(ctrl.input_change, "['csr_bins']"),
                        )
