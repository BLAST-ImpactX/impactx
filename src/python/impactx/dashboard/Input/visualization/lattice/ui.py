"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import setup_server, vuetify
from ... import CardBase, CardComponents
server, state, ctrl = setup_server()

ELEMENT_COLOR_MAP = {
    "drift": "blue lighten-2",
    "quad": "red darken-1",
    "monitor": "grey darken-2",
}

def get_element_color(name: str) -> str:
    """
    Determine Vuetify color for an element based on its name.
    """
    clean_name = name.lower()
    for element_key, color in ELEMENT_COLOR_MAP.items():
        if element_key in clean_name:
            return color
    return "grey lighten-1"


@state.change("selected_lattice_list")
def on_lattice_list_change(**kwargs):
    for element in state.selected_lattice_list:
        element["color"] = get_element_color(element["name"])

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
                            color=("element.color",),
                            text_color="white",
                            elevation=2,
                        ):
                            pass