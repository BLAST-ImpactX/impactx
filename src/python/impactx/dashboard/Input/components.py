from typing import Optional

from .. import html, setup_server, vuetify
from .defaults import TooltipDefaults
from .generalFunctions import generalFunctions

server, state, ctrl = setup_server()

state.documentation_drawer_open = False
state.documentation_url = ""


class CardComponents:
    """
    Class contains staticmethods to build
    card components using Vuetify's VCard.
    """

    @staticmethod
    def input_header(section_name: str, additional_components=None) -> None:
        """
        Creates a standardized header look for inputs.

        :param section_name: The name for the input section.
        """

        documentation_name = section_name.lower().replace(" ", "_")
        with vuetify.VCardTitle(section_name):
            vuetify.VSpacer()
            if additional_components:
                additional_components()
            CardComponents.refresh_icon(documentation_name)
            CardComponents.documentation_icon(documentation_name)
        vuetify.VDivider()

    @staticmethod
    def documentation_icon(section_name: str) -> vuetify.VIcon:
        """
        Takes user to input section's documentation.

        :param section_name: The name for the input section.
        """

        return vuetify.VIcon(
            "mdi-information",
            style="color: #00313C;",
            click=lambda: generalFunctions.open_documentation(section_name),
        )

    @staticmethod
    def refresh_icon(section_name: str) -> vuetify.VIcon:
        """
        Resets input values to default.

        :param section_name: The name for the input section.
        """

        return vuetify.VIcon(
            "mdi-refresh",
            style="color: #00313C;",
            click=lambda: generalFunctions.reset_inputs(section_name),
        )


class InputComponents:
    """
    Class contains staticmethod to create
    input-related Vuetify components.
    """

    DENSE = True

    @staticmethod
    def select(
        label: str,
        v_model_name: Optional[str] = None,
        items: Optional[list] = None,
        **kwargs,
    ) -> vuetify.VSelect:
        """
        Creates a Vuetify VSelect component with
        pre-filled components.

        :param label: Display label
        :param v_model_name: v_model binding name. Optional, as default name
        created otherwise with label name.
        :param items: Items list override
        """

        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        if items is None:
            items = (
                generalFunctions.get_default(f"{v_model_name}_list", "default_values"),
            )

        with vuetify.VTooltip(**TooltipDefaults.TOOLTIP_STYLE):
            with vuetify.Template(v_slot_activator="{ on, attrs }"):
                vuetify.VSelect(
                    label=label,
                    v_model=(v_model_name,),
                    items=items,
                    dense=True,
                    **kwargs,
                    v_on="on",
                    v_bind="attrs",
                )
            html.Span(TooltipDefaults.TOOLTIP.get(v_model_name))

    @staticmethod
    def text_field(
        label: str, v_model_name: Optional[str] = None, **kwargs
    ) -> vuetify.VTextField:
        """
        Creates a Vuetify VTextField component with the following default components:
        - error_message state: It's init value is an empty list.
        - step: The step value of the input (either set in defaults.py), or
          by default is set to 1.
        - suffix: The unit of the input (either set in defauts.py), or
          by default is empty.
        - type: set to 'number' to only allow a numeric input.
        - dense: set to 'true' to minimize space usage.

        :param label: Display label
        :param v_model_name: v_model binding name. Optional, as default name
        created otherwise with label name.
        """

        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        with vuetify.VTooltip(**TooltipDefaults.TOOLTIP_STYLE):
            with vuetify.Template(v_slot_activator="{ on, attrs }"):
                vuetify.VTextField(
                    label=label,
                    v_model=(v_model_name,),
                    error_messages=(f"{v_model_name}_error_message", []),
                    type="number",
                    step=generalFunctions.get_default(f"{v_model_name}", "steps"),
                    suffix=generalFunctions.get_default(f"{v_model_name}", "units"),
                    __properties=["step"],
                    dense=True,
                    v_on="on",
                    v_bind="attrs",
                    **kwargs,
                )
            html.Span(TooltipDefaults.TOOLTIP.get(v_model_name))


class NavigationComponents:
    """
    Class contains staticmethods to create
    navigation-related Vuetify components.
    """

    @staticmethod
    def create_route(route_title: str, mdi_icon: str) -> None:
        """
        Creates a route with specified title and icon.

        :param route_title: The title for the route
        :param mdi_icon: The MDI icon name to display
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
    def create_dialog_tabs(name: str, num_tabs: int, tab_names: list[str]) -> None:
        """
        Creates a tabbed dialog interface.

        :param name: The base name for the tab group
        :param num_tabs: Number of tabs to create
        :param tab_names: List of names for each tab
        """
        if len(tab_names) != num_tabs:
            raise ValueError("Number of tab names must match number of tabs_names")

        with vuetify.VCard():
            with vuetify.VTabs(v_model=(f"{name}", 0)):
                for tab_name in tab_names:
                    vuetify.VTab(tab_name)
            vuetify.VDivider()

    @staticmethod
    def create_documentation_drawer():
        with vuetify.VNavigationDrawer(
            v_model=("documentation_drawer_open",),
            absolute=True,
            right=True,
            hide_overlay=True,
            style="width: 30vw; top: 64px !important;  position: fixed;",
        ):
            with vuetify.VContainer(
                fluid=True,
                classes="pa-0 fill-height",
            ):
                html.Iframe(
                    src=("documentation_url",),
                    style="width: 100%; height: 100%; border: none;",
                )
