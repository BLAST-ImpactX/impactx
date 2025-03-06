"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardBase, CardComponents, InputComponents

server, state, ctrl = setup_server()


class csrConfiguration(CardBase):
    HEADER_NAME = "CSR"

    def __init__(self):
        super().__init__()

    def card_content(self):
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(self.HEADER_NAME)
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
