from typing import Optional

from . import setup_server, vuetify
from .generalFunctions import generalFunctions

server, state, ctrl = setup_server()


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
            click=lambda: generalFunctions.documentation(section_name),
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

        return vuetify.VSelect(
            label=label,
            v_model=(v_model_name,),
            items=items,
            dense=True,
            **kwargs,
        )

    @staticmethod
    def text_field(
        label: str, v_model_name: Optional[str] = None, **kwargs
    ) -> vuetify.VTextField:
        """
        Creates a Vuetify VTextField component with
        pre-filled components.

        :param label: Display label
        :param v_model_name: v_model binding name. Optional, as default name
        created otherwise with label name.
        """

        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        return vuetify.VTextField(
            label=label,
            v_model=(v_model_name,),
            error_messages=(f"{v_model_name}_error_message", []),
            type="number",
            step=generalFunctions.get_default(f"{v_model_name}", "steps"),
            suffix=generalFunctions.get_default(f"{v_model_name}", "units"),
            __properties=["step"],
            dense=True,
            **kwargs,
        )
