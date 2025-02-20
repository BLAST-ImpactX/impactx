from ... import setup_server, vuetify

server, state, ctrl = setup_server()


class LatticeConfigurationHelper:
    """
    Helper class to build the Lattice Configuration section of the dashboard
    """

    def expand_configuration() -> vuetify.VBtn:
        """
        A button which expands/closes the lattice configuration.
        """

        with vuetify.VBtn(
            color="primary",
            click="expand_configuration = !expand_configuration",
            icon=True,
            small=True,
        ):
            vuetify.VIcon(
                v_text=("expand_configuration ? 'mdi-close' : 'mdi-arrow-expand'",),
            )

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
