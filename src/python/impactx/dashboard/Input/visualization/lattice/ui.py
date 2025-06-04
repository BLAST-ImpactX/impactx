"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import html, setup_server, vuetify
from ... import CardBase, NavigationComponents
from . import Dialogs, StatComponents, StatUtils

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

def _update_statistics():
    """
    Update statistics based on the current selected lattice elements.
    """
    for element in state.selected_lattice_list:
        element["color"] = get_element_color(element["name"])

    state.total_elements = len(state.selected_lattice_list)
    state.periods = 1 if state.total_elements > 0 else 0
    state.total_steps = StatUtils.update_total_steps()
    state.element_counts = StatUtils.update_element_counts()
    StatUtils.update_length_statistics()

@state.change("selected_lattice_list")
def on_lattice_list_change(**kwargs):
    _update_statistics()

class LatticeVisualizer(CardBase):

    def __init__(self):
        super().__init__()
        
    def card_content(self):
        with vuetify.VDialog(
            v_model=("lattice_visualizer_dialog_settings", False), max_width="33.33vw"
        ):
            self.dialog_settings()

        with vuetify.VCard():
            with vuetify.VCard(
                classes="d-flex flex-column",
                style="min-height: 3.75rem; margin-bottom: 20px;",
                color="grey lighten-4",
                elevation=2,
            ):
                # create custom header over using component in CardComponents
                with vuetify.VCardTitle(classes="d-flex align-center"):
                    html.Div("Lattice Statistics")
                    vuetify.VSpacer()
                    StatComponents.settings()
                StatComponents.statistics()

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

    @staticmethod
    def dialog_settings():
        dialog_name = "lattice_visualizer_dialog_tab_settings"

        with NavigationComponents.create_dialog_tabs(
            dialog_name, 2, ["Element Colors", "General Settings"]
        ):
            with vuetify.VTabsWindow(v_model=(dialog_name, 0)):
                with vuetify.VTabsWindowItem():
                    Dialogs.element_colors_tab()
                with vuetify.VTabsWindowItem():
                    Dialogs.general_settings_tab()
