"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import vuetify

from ..trame_setup import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


class TrameFunctions:
    """
    Contains functions containing Vuetify
    components.
    """

    @staticmethod
    def create_route(route_title, mdi_icon):
        """
        Creates a route with a specified title and icon.
        :param route_title: The title of the route.
        :param mdi_icon: The icon to be used for the route.
        """

        state[route_title] = False  # Does not display route by default

        to = f"/{route_title}"
        click = f"{route_title} = true"

        with vuetify.VListItem(to=to, click=click):
            with vuetify.VListItemIcon():
                vuetify.VIcon(mdi_icon)
            with vuetify.VListItemContent():
                vuetify.VListItemTitle(route_title)

    @staticmethod
    def create_dialog_tabs(name: str, num_tabs: int, tab_names: list[str]):
        if len(tab_names) != num_tabs:
            raise ValueError("Number of tab names must match number of tabs_names")

        with vuetify.VCard():
            with vuetify.VTabs(v_model=(f"{name}", 0)):
                for tab_name in tab_names:
                    vuetify.VTab(tab_name)
            vuetify.VDivider()
