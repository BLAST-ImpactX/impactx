from ... import html, setup_server, vuetify
from ..defaults import UIDefaults
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

        self.card_props = {"elevation": 2, "style": self.collapsable}

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

        with vuetify.VCardTitle(
            section_name,
            classes="d-flex align-center flex-wrap",
            style="min-height: 3.75rem;",
        ):
            vuetify.VSpacer()
            with html.Div(classes="d-flex", gap="2px"):
                render_components("start")
                CardComponents.documentation_button(section_name_cleaned)
                CardComponents.refresh_button(section_name_cleaned)
                CardComponents.collapse_button(section_name_cleaned)
                CardComponents.expand_button(section_name_cleaned)
                render_components("end")
        vuetify.VDivider()

    @staticmethod
    def card_button(
        icon_name, color="primary", dynamic_condition=None, description=None, **kwargs
    ) -> vuetify.VBtn:
        """
        Create a Vuetify VBtn containing an icon.

        :param icon_name: A string for a static icon, or a list/tuple of two strings for conditional rendering.
        :param color: The button color.
        :param dynamic_condition: A Vue expression that determines which icon to display when `icon_name` is a list/tuple.
        :param kwargs: Extra keyword arguments for the VBtn component.
        """

        with vuetify.VTooltip(location="bottom", text=description):
            with vuetify.Template(v_slot_activator="{ props }"):
                with vuetify.VBtn(
                    color=color,
                    icon=True,
                    density="compact",
                    variant="text",
                    v_bind="props",
                    **kwargs,
                ):
                    if isinstance(icon_name, (list, tuple)):
                        with vuetify.Template(v_if=dynamic_condition):
                            vuetify.VIcon(icon_name[1])
                        with vuetify.Template(v_else=True):
                            vuetify.VIcon(icon_name[0])
                    else:
                        vuetify.VIcon(icon_name)

    @staticmethod
    def documentation_button(section_name: str) -> vuetify.VBtn:
        """
        Takes user to input section's documentation.

        :param section_name: The name for the input section.
        """

        CardComponents.card_button(
            "mdi-information",
            color="#00313C",
            click=lambda: generalFunctions.open_documentation(section_name),
            description="Documentation",
        )

    @staticmethod
    def refresh_button(section_name: str) -> vuetify.VBtn:
        """
        Resets input values to default.

        :param section_name: The name for the input section.
        """

        CardComponents.card_button(
            "mdi-refresh",
            color="#00313C",
            click=lambda: generalFunctions.reset_inputs(section_name),
            description="Reset",
        )

    @staticmethod
    def expand_button(section_name: str) -> vuetify.VBtn:
        """
        A button which expands/closes the given card configuration.

        :param section_name: The name for the input section.
        """

        expand_state = f"expand_{section_name}"

        CardComponents.card_button(
            ["mdi-arrow-expand", "mdi-close"],
            click=f"{expand_state} = !{expand_state}",
            dynamic_condition=expand_state,
            description="Expand",
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

        CardComponents.card_button(
            ["mdi-chevron-up", "mdi-chevron-down"],
            click=f"{collapsed_state_name} = !{collapsed_state_name}",
            dynamic_condition=collapsed_state_name,
            description="Collapse",
        )
