"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardComponents, InputComponents, UIDefaults

server, state, ctrl = setup_server()


class csrConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for CSR.
        """

        with vuetify.VCard(**UIDefaults.card_sizing):
            CardComponents.input_header("CSR")
            with vuetify.VCardText(**UIDefaults.card_text_overflow):
                with vuetify.VRow(**UIDefaults.row_style):
                    with vuetify.VCol():
                        InputComponents.select(
                            label="Particle Shape",
                        )
                with vuetify.VRow(**UIDefaults.row_style):
                    with vuetify.VCol():
                        InputComponents.text_field(
                            label="CSR Bins",
                        )
