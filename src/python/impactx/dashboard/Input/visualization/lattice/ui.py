"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import setup_server, vuetify
from ... import CardBase, CardComponents
server, state, ctrl = setup_server()


@state.change("selected_lattice_list")
def on_lattice_list_change(**kwargs):
    print(f"lattice list changed")

class LatticeVisualizer(CardBase):
    HEADER_NAME = "Lattice Visualizer"

    def __init__(self):
        super().__init__()

    def card_content(self):
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(self.HEADER_NAME)
            with vuetify.VCardText():
                with vuetify.VRow():
                    with vuetify.VCol(
                        v_for="(element, index) in selected_lattice_list",
                    ):
                        with vuetify.VCard(
                            text=("element.name",),
                            class_="rounded-xl pa-2 text-center",
                            color="primary",
                            text_color="white",
                            elevation=2,
                        ):
                            pass