"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from trame.widgets import vuetify

from ..Input.generalFunctions import generalFunctions
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
    def select(label, v_model_name=None, items=None, **kwargs):
        # in place for now as some variables are not in same format
        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        if items is None:
            items = (
                generalFunctions.get_default(f"{v_model_name}_list", "default_values"),
            )

        vuetify.VSelect(
            label=label,
            v_model=(f"{v_model_name}",),
            items=items,
            dense=True,
            **kwargs,
        )

    @staticmethod
    def text_field(label, v_model_name=None, **kwargs):
        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        vuetify.VTextField(
            label=label,
            v_model=(f"{v_model_name}",),
            error_messages=(f"{v_model_name}_error_message", []),
            type="number",
            step=generalFunctions.get_default(f"{v_model_name}", "steps"),
            suffix=generalFunctions.get_default(f"{v_model_name}", "units"),
            __properties=["step"],
            dense=True,
            **kwargs,
        )

    @staticmethod
    def documentation_icon(section_name, **kwargs):
        vuetify.VIcon(
            "mdi-information",
            style="color: #00313C;",
            click=lambda: generalFunctions.documentation(section_name),
        )

    @staticmethod
    def refresh_icon(section_name):
        """
        Creates a standardized refresh button.
        :param reset_function: The reset function to call when clicked.
        """
        vuetify.VIcon(
            "mdi-refresh",
            style="color: #00313C;",
            click=lambda: generalFunctions.reset_inputs(section_name),
        )

    @staticmethod
    def input_section_header(section_name, additional_components=None):
        documentation_name = section_name.lower().replace(" ", "_")
        with vuetify.VCardTitle(section_name):
            vuetify.VSpacer()
            if additional_components:
                additional_components()
            TrameFunctions.refresh_icon(documentation_name)
            TrameFunctions.documentation_icon(documentation_name)
        vuetify.VDivider()

    @staticmethod
    def create_dialog_tabs(name: str, num_tabs: int, tab_names: list[str]):
        if len(tab_names) != num_tabs:
            raise ValueError("Number of tab names must match number of tabs_names")

        with vuetify.VCard():
            with vuetify.VTabs(v_model=(f"{name}", 0)):
                for tab_name in tab_names:
                    vuetify.VTab(tab_name)
            vuetify.VDivider()
