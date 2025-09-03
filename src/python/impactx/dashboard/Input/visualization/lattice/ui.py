"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import html, state, vuetify
from ....Input.components import CardBase
from . import StatComponents, StatUtils


def _update_statistics() -> None:
    """
    Update statistics based on the current selected lattice elements.
    """
    state.total_elements = len(state.selected_lattice_list)
    state.total_steps = StatUtils.update_total_steps()
    state.element_counts = StatUtils.update_element_counts()
    StatUtils.update_length_statistics()


@state.change("selected_lattice_list")
def on_lattice_list_change(**kwargs):
    _update_statistics()


class LatticeVisualizer(CardBase):
    """
    Displays the lattice visualizer section on the inputs page of the dashboard.
    """

    HEADER_NAME = "Lattice Statistics"
    SUPPRESS_DOC_WARNING = True

    def __init__(self):
        super().__init__()

    def card_content(self):
        """
        The content of the lattice visualizer.
        """

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
                StatComponents.statistics()
