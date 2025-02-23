from typing import Optional

from ... import html, setup_server, vuetify
from ..defaults import TooltipDefaults, UIDefaults
from ..generalFunctions import generalFunctions

server, state, ctrl = setup_server()

state.documentation_drawer_open = False
state.documentation_url = ""


def clean_name(section_name):
    return section_name.lower().replace(" ", "_")


class CardBase(UIDefaults):
    HEADER_NAME = "Base Section"

    def __init__(self):
        self.header = self.HEADER_NAME.lower().replace(" ", "_")
        self.collapsable = (f"collapse_{self.header}_height",)

    def card(self):
        """
        Creates UI content for a section.
        """

        self.init_dialog(self.HEADER_NAME, self.card_content)
        self.card_content()

    def card_content(self):
        raise NotImplementedError("Card must contain card_content.")

    @staticmethod
    def init_dialog(section_name: str, content_callback) -> None:
        """
        Renders the expansion dialog UI for the input sections card.
        Only runs once, when the section's card is built.
        """

        section_name_cleaned = clean_name(section_name)
        expand_state_name = f"expand_{section_name_cleaned}"

        setattr(state, expand_state_name, False)

        with vuetify.VDialog(v_model=(expand_state_name,), width="fit-content"):
            with vuetify.VCard():
                content_callback()


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

        section_name_cleaned = clean_name(section_name)

        def render_components(position: str):
            if additional_components and position in additional_components:
                additional_components[position]()

        with vuetify.VCardTitle(section_name):
            vuetify.VSpacer()
            render_components("start")
            CardComponents.refresh_icon(section_name_cleaned)
            CardComponents.documentation_icon(section_name_cleaned)
            CardComponents.collapse_button(section_name_cleaned)
            CardComponents.expand_button(section_name_cleaned)
            render_components("end")
        vuetify.VDivider()

    @staticmethod
    def documentation_icon(section_name: str) -> vuetify.VBtn:
        """
        Takes user to input section's documentation.

        :param section_name: The name for the input section.
        """

        with vuetify.VBtn(
            style="color: #00313C;",
            click=lambda: generalFunctions.open_documentation(section_name),
            icon=True,
            small=True,
        ):
            vuetify.VIcon(
                "mdi-information",
            )

    @staticmethod
    def refresh_icon(section_name: str) -> vuetify.VBtn:
        """
        Resets input values to default.

        :param section_name: The name for the input section.
        """

        with vuetify.VBtn(
            style="color: #00313C;",
            click=lambda: generalFunctions.reset_inputs(section_name),
            icon=True,
            small=True,
        ):
            vuetify.VIcon("mdi-refresh")

    @staticmethod
    def expand_button(section_name: str) -> vuetify.VBtn:
        """
        A button which expands/closes the given card configuration.

        :param section_name: The name for the input section.
        """

        with vuetify.VBtn(
            color="primary",
            click=f"expand_{section_name} = !expand_{section_name}",
            icon=True,
            small=True,
        ):
            vuetify.VIcon(
                v_text=(f"expand_{section_name} ? 'mdi-close' : 'mdi-arrow-expand'",)
            )

    @staticmethod
    def collapse_button(section_name: str) -> vuetify.VBtn:
        """
        A button which collapses the given cards inputs.

        :param section_name: The name for the input section.
        """
        section_name_cleaned = clean_name(section_name)
        collapsed_state_name = f"collapse_{section_name_cleaned}"

        setattr(state, collapsed_state_name, False)

        with vuetify.VBtn(
            color="primary",
            click=f"collapse_{section_name_cleaned} = !collapse_{section_name_cleaned}",
            icon=True,
            small=True,
        ):
            vuetify.VIcon(
                v_text=(
                    f"collapse_{section_name_cleaned} ? 'mdi-chevron-down' : 'mdi-chevron-up'",
                )
            )
