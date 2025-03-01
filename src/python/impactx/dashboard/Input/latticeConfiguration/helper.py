from ... import setup_server, vuetify
from .. import CardComponents

server, state, ctrl = setup_server()


class LatticeConfigurationHelper:
    """
    Helper class to build the Lattice Configuration section of the dashboard
    """

    BUTTON_COLOR = "grey-darken-2"
    BUTTON_COLOR_LIGHTER = "grey-darken-1"

    @staticmethod
    def settings() -> vuetify.VBtn:
        """
        A button which opens the lattice configuration settings.
        """

        CardComponents.card_button(
            "mdi-cog",
            color=LatticeConfigurationHelper.BUTTON_COLOR,
            click="lattice_configuration_dialog_settings = true",
            documentation="Settings",
        )

    @staticmethod
    def move_element_up() -> vuetify.VBtn:
        """
        A button which allows the dashboard user to
        move a lattice element's index upward.
        """

        CardComponents.card_button(
            "mdi-menu-up",
            color=LatticeConfigurationHelper.BUTTON_COLOR_LIGHTER,
            click=(ctrl.move_latticeElementIndex_up, "[index]"),
        )

    @staticmethod
    def move_element_down() -> vuetify.VBtn:
        """
        A button which allows the dashboard user to
        move a lattice element's index downward.
        """

        CardComponents.card_button(
            "mdi-menu-down",
            color=LatticeConfigurationHelper.BUTTON_COLOR_LIGHTER,
            click=(ctrl.move_latticeElementIndex_down, "[index]"),
        )

    @staticmethod
    def delete_element() -> vuetify.VBtn:
        """
        A button which allows the dashboard user to
        move a lattice element's index downward.
        """

        CardComponents.card_button(
            "mdi-delete",
            color=LatticeConfigurationHelper.BUTTON_COLOR,
            click=(ctrl.deleteLatticeElement, "[index]"),
        )
