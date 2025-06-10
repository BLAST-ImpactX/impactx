"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import html, setup_server, vuetify
from trame.widgets import plotly

from ... import CardBase, NavigationComponents
from . import Dialogs, StatComponents, StatUtils
from .visualization.plot import lattice_visualizer

server, state, ctrl = setup_server()


def _update_statistics() -> None:
    """
    Update statistics based on the current selected lattice elements.
    """
    state.total_elements = len(state.selected_lattice_list)
    state.periods = 1 if state.total_elements > 0 else 0
    state.total_steps = StatUtils.update_total_steps()
    state.element_counts = StatUtils.update_element_counts()
    StatUtils.update_length_statistics()

def _update_lattice_visualization() -> None:
    """
    Updates the plotly figure with an updated lattice visualization.
    """
    ctrl.lattice_figure_update(lattice_visualizer())

@state.change("selected_lattice_list")
def on_lattice_list_change(**kwargs):
    _update_statistics()
    _update_lattice_visualization()

class LatticeVisualizer(CardBase):
    """
    Displays the lattice visualizer section on the inputs page of the dashboard.
    """

    def __init__(self):
        super().__init__()
        
    def card_content(self):
        """
        The content of the lattice visualizer.
        """
        with vuetify.VDialog(
            v_model=("lattice_visualizer_dialog_settings", False), max_width="33.33vw"
        ):
            self.dialog_settings()

        with vuetify.VCard():
            with vuetify.VCard(
                classes="d-flex flex-column",
                style="min-height: 3.75rem; margin-bottom: 20px;",
                color="#002949",
                elevation=2,
            ):
                # create custom header over using component in CardComponents
                with vuetify.VCardTitle(classes="d-flex align-center"):
                    html.Div("Lattice Statistics")
                    vuetify.VSpacer()
                    Dialogs.settings()
                StatComponents.statistics()
            
            with vuetify.VCard(color="#002949"):
                with vuetify.VCardText():
                    ctrl.lattice_figure_update = plotly.Figure(
                        display_mode_bar="true",
                        style="width: 100%; height: 50vh"
                    ).update

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
