"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import ctrl, html, state, vuetify
from ..utils import DEFAULT_HEADERS, UNSELECTABLE_HEADERS
from .components import components

state.over_s_header_search = ""
state.available_headers = []
state.available_headers_on_ui = []
state.unselectable_headers = UNSELECTABLE_HEADERS
state.no_results = False
state.no_available_items = False


def _filter_headers_by_query() -> list[str]:
    """
    Filters the available headers based on the search query.
    """

    headers = state.available_headers
    query = state.over_s_header_search

    query = (query or "").lower()
    if query:
        return [header for header in headers if query in header.lower()]
    state.dirty("available_headers_on_ui")
    return headers


def _sync_headers():
    """
    Ensures that the selected headers are in sync with the available headers.

    The available headers are those that are selectable
    and not already selected by the user.
    """

    selectable_keys = [item["key"] for item in state.selectable_headers or [] if item]

    state.available_headers = [
        h for h in selectable_keys if h not in state.selected_headers
    ]
    state.available_headers_on_ui = _filter_headers_by_query()
    state.dirty("available_headers")


@state.change("selected_headers", "selectable_headers")
def on_header_state_change(**kwargs):
    """
    Called when the set of selectable headers is modified
    by the user or after the initial simulation run
    for the 'Plot Over S' visualization.
    """
    _sync_headers()


@ctrl.add("reset_selected_headers")
def reset_selected_headers():
    """
    Resets the selected y-axis headers to the default state.
    """
    state.selected_headers = list(DEFAULT_HEADERS)


@ctrl.add("add_selected_header")
def add_selected_header(item):
    state.selected_headers = state.selected_headers + [item]


@ctrl.add("remove_selected_header")
def remove_selected_header(item):
    state.selected_headers = [h for h in state.selected_headers if h != item]


@state.change("over_s_header_search")
def on_over_s_header_search_change(**kwargs):
    state.available_headers_on_ui = _filter_headers_by_query()


@state.change("available_headers_on_ui", "over_s_header_search")
def update_empty_states(**kwargs):
    has_search = bool(state.over_s_header_search)
    state.no_results = has_search and not state.available_headers_on_ui
    state.no_available_items = not has_search and not state.available_headers_on_ui
    state.dirty("no_results", "no_available_items")


class OverSHeaderSelector:
    def selector(self):
        """
        Displays the selector for header selection.
        """

        with vuetify.VCard(
            elevation=2,
            rounded="lg",
            classes="d-flex flex-column",
            style="height: 87vh;",
        ):
            with vuetify.VCardTitle(classes="pa- d-flex"):
                html.Span(
                    "Select headers to plot", classes="text-subtitle-1 font-weight-bold"
                )
                vuetify.VSpacer()
                components.reset()
            vuetify.VDivider()

            with html.Div(classes="flex-grow-1 overflow-hidden"):
                with vuetify.VRow(no_gutters=True, classes="fill-height"):
                    components.selected_headers_column()
                    vuetify.VDivider(vertical=True)
                    components.available_headers_column()


over_s_selector = OverSHeaderSelector()
