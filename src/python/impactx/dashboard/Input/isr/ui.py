"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardBase, CardComponents, InputComponents

server, state, ctrl = setup_server()


class isrConfiguration(CardBase):
    HEADER_NAME = "ISR"

    def __init__(self):
        super().__init__()

    def card_content(self):
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(self.HEADER_NAME)
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol():
                        InputComponents.select(
                            label="ISR Order",
                        )
