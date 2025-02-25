from ... import setup_server, vuetify

server, state, ctrl = setup_server()


class LatticeConfigurationHelper:
    """
    Helper class to build the Lattice Configuration section of the dashboard
    """

    @staticmethod
    def settings() -> vuetify.VBtn:
        """
        A button which opens the lattice configuration settings.
        """

        with vuetify.VBtn(
            click="lattice_configuration_dialog_settings = true",
            icon=True,
            small=True,
        ):
            vuetify.VIcon("mdi-cog")

    @staticmethod
    def move_element_up() -> vuetify.VBtn:
        """
        A button which allows the dashboard user to
        move a lattice element's index upward.
        """

        with vuetify.VBtn(
            click=(ctrl.move_latticeElementIndex_up, "[index]"),
            icon=True,
            small=True,
        ):
            vuetify.VIcon(
                "mdi-menu-up",
            )

    @staticmethod
    def move_element_down() -> vuetify.VBtn:
        """
        A button which allows the dashboard user to
        move a lattice element's index downward.
        """

        with vuetify.VBtn(
            click=(ctrl.move_latticeElementIndex_down, "[index]"),
            icon=True,
            small=True,
        ):
            vuetify.VIcon(
                "mdi-menu-down",
            )

    @staticmethod
    def delete_element() -> vuetify.VBtn:
        """
        A button which allows the dashboard user to
        move a lattice element's index downward.
        """

        with vuetify.VBtn(
            click=(ctrl.deleteLatticeElement, "[index]"),
            icon=True,
            small=True,
        ):
            vuetify.VIcon(
                "mdi-delete",
            )
