"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardComponents, InputComponents

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
                        )
