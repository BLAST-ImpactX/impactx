from ... import html, setup_server, vuetify

server, state, ctrl = setup_server()

state.documentation_drawer_open = False
state.documentation_url = ""


class NavigationComponents:
    """
    Class contains staticmethods to create
    navigation-related Vuetify components.
    """

    @staticmethod
    def create_route(route_title: str, mdi_icon: str) -> vuetify.VListItem:
        """
        Creates a route with specified title and icon.

        :param route_title: The title for the route
        :param mdi_icon: The MDI icon name to display
        """
        state[route_title] = False  # Does not display route by default

        to = f"/{route_title}"
        click = f"{route_title} = true"

        return vuetify.VListItem(
            to=to,
            click=click,
            prepend_icon=mdi_icon,
            title=route_title,
            style="height: 3rem;",
        )

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

        card = vuetify.VCard()
        with card:
            with vuetify.VTabs(v_model=(f"{name}", 0)):
                for tab_name in tab_names:
                    vuetify.VTab(tab_name)
            vuetify.VDivider()
        return card

    @staticmethod
    def create_documentation_drawer():
        with vuetify.VNavigationDrawer(
            model_value=("documentation_drawer_open",),
            location="right",
            temporary=True,
            scrim=False,
            style="z-index: 10",
            width=500,
        ):
            with vuetify.VContainer(
                classes="pa-0 fill-height",
            ):
                html.Iframe(
                    src=("documentation_url",),
                    style="width: 100%; height: 100%; border: none;",
                )
