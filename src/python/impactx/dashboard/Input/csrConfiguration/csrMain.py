"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardComponents, InputComponents, UIDefaults

server, state, ctrl = setup_server()


class csrConfiguration(UIDefaults):
    def __init__(self):
        super().__init__()

    def card(self):
        """
        Creates UI content for CSR.
        """

        with vuetify.VCard(style=("card_style",)):
            CardComponents.input_header("CSR")
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol():
                        InputComponents.select(
                            label="Particle Shape",
                        )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol():
                        InputComponents.text_field(
                            label="CSR Bins",
                        )
