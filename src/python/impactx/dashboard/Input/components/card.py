from typing import List, Optional, Union

from ... import html, state, vuetify
from ..defaults import DashboardDefaults, UIDefaults
from ..utils import GeneralFunctions

state.documentation_drawer_open = False
state.documentation_url = ""

_missing_docs = set()


def clean_name(section_name):
    return GeneralFunctions.normalize_for_v_model(section_name)


class CardBase(UIDefaults):
    HEADER_NAME = "Base Section"

    def __init__(self):
        self.header = GeneralFunctions.normalize_for_v_model(self.HEADER_NAME)
        self.collapsable = (f"collapse_{self.header}_height",)

        self.card_props = {"elevation": 2, "style": self.collapsable}

    def card(self):
        """
        Creates UI content for a section.
        """
        # Allow subclasses to suppress missing documentation warning
        if not getattr(self, "SUPPRESS_DOC_WARNING", False):
            if (
                self.header not in DashboardDefaults.DOCUMENTATION
                and self.header not in _missing_docs
            ):
                print(
                    f"WARNING: Card '{self.header}' has no doc link in DashboardDefaults.DOCUMENTATION"
                )
                _missing_docs.add(self.header)

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
        icon_name: Union[str, List[str]],
        color: str = "primary",
        dynamic_condition: Optional[str] = None,
        description: Optional[Union[str, List[str]]] = None,
        density: str = "compact",
        variant: str = "text",
        **kwargs,
    ) -> vuetify.VBtn:
        """
        Creates a Vuetify VBtn as an icon button. Can be dynamically toggled
        between two states using 'dynamic_condition'.

        :param icon_name: A string or a list of two strings for conditional rendering of the button icon.
        :param color: The button color.
        :param dynamic_condition: A Vue state variable (boolean) that controls which value to show in 'icon_name' and 'description'.
        :param description: A string or a list of two strings for conditional tooltip text.
        :param kwargs: Extra keyword arguments for the component.
        """

        def validate_dynamic_condition(prop_value: List[str], prop_name: str) -> None:
            """
            Ensure dynamic_condition components are a list of exactly 2 strings for dynamic toggling (e.g., expand/collapse).
            """

            if not isinstance(prop_value, (list, tuple)):
                raise ValueError(
                    f"When dynamic_condition is set, {prop_name} must be a list of exactly 2 strings"
                )
            if len(prop_value) != 2:
                raise ValueError(
                    f"When dynamic_condition is set, {prop_name} must contain exactly 2 elements"
                )
            if not all(isinstance(item, str) for item in prop_value):
                raise ValueError(
                    f"When dynamic_condition is set, all elements in {prop_name} must be strings"
                )

        if dynamic_condition:
            if description is None:
                raise ValueError(
                    "When dynamic_condition is set, 'description' must be provided and cannot be None."
                )
            validate_dynamic_condition(icon_name, "icon_name")
            validate_dynamic_condition(description, "description")

        if dynamic_condition:
            tooltip_text = (
                f"{dynamic_condition} ? '{description[1]}' : '{description[0]}'",
            )
        else:
            tooltip_text = description

        with vuetify.VTooltip(location="bottom", text=tooltip_text):
            with vuetify.Template(v_slot_activator="{ props }"):
                with vuetify.VBtn(
                    color=color,
                    icon=True,
                    density=density,
                    variant=variant,
                    v_bind="props",
                    **kwargs,
                ):
                    if dynamic_condition:
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
            click=lambda: GeneralFunctions.open_documentation(section_name),
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
            click=lambda: GeneralFunctions.reset_inputs(section_name),
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
            description=["Expand", "Close"],
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
            description=["Minimize", "Show"],
        )
