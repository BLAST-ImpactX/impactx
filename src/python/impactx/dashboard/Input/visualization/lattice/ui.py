"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import setup_server, vuetify
from ... import CardBase, CardComponents, NavigationComponents
from . import LatticeVisualizerComponents, LatticeVisualizerDialogs, LatticeVisualizerUtils

server, state, ctrl = setup_server()

utils = LatticeVisualizerUtils

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

    state.total_elements = len(state.selected_lattice_list)
    state.total_steps = utils.update_total_steps()
    state.element_counts = utils.update_element_counts()
    utils.update_length_statistics()

class LatticeVisualizer(CardBase):
    HEADER_NAME = "Lattice Visualizer"
    components = LatticeVisualizerComponents

    def __init__(self):
        super().__init__()
        
    def card_content(self):
        with vuetify.VDialog(
            v_model=("lattice_visualizer_dialog_settings", False), max_width="33.33vw"
        ):
            self.dialog_settings()

        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(
                self.HEADER_NAME,
                additional_components={"end": self.components.settings},
            )
            with vuetify.VCardText(): 
                self.components.statistics()
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

    @staticmethod
    def dialog_settings():
        dialog_name = "lattice_visualizer_dialog_tab_settings"

        with NavigationComponents.create_dialog_tabs(
            dialog_name, 2, ["Element Colors", "General Settings"]
        ):
            with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
                with vuetify.VTabsWindowItem():
                    LatticeVisualizerDialogs.element_colors_tab()
                with vuetify.VTabsWindowItem():
                    LatticeVisualizerDialogs.general_settings_tab()
