from typing import Optional

from ... import html, setup_server, vuetify
from ..defaults import TooltipDefaults
from ..generalFunctions import generalFunctions

server, state, ctrl = setup_server()


class InputComponents:
    """
    Class contains staticmethod to create
    input-related Vuetify components.
    """

    @staticmethod
    def _build_component(
        vuetify_component,
        label: str,
        v_model_name: Optional[str] = None,
        **component_kwargs,
    ) -> None:
        """
        Helper to build a component with common properties and tooltip/template wrappers.
        """
        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        if "items" in component_kwargs and component_kwargs["items"] is None:
            component_kwargs["items"] = (
                generalFunctions.get_default(f"{v_model_name}_list", "default_values"),
            )

        common_props = {
            "label": label,
            "v_model": (v_model_name,),
            "dense": True,
            "v_on": "on",
            "v_bind": "attrs",
        }
        props = {**common_props, **component_kwargs}

        with vuetify.VTooltip(**TooltipDefaults.TOOLTIP_STYLE):
            with vuetify.Template(v_slot_activator="{ on, attrs }"):
                vuetify_component(**props)
            html.Span(TooltipDefaults.TOOLTIP.get(v_model_name))

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
        InputComponents._build_component(
            vuetify.VSelect, label, v_model_name, items=items, **kwargs
        )

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
        computed_v_model = (
            v_model_name
            if v_model_name is not None
            else label.lower().replace(" ", "_")
        )
        InputComponents._build_component(
            vuetify.VTextField,
            label,
            v_model_name,
            error_messages=(f"{computed_v_model}_error_message", []),
            type="number",
            step=generalFunctions.get_default(computed_v_model, "steps"),
            suffix=generalFunctions.get_default(computed_v_model, "units"),
            __properties=["step"],
            **kwargs,
        )

    @staticmethod
    def checkbox(
        label: str, v_model_name: Optional[str] = None, **kwargs
    ) -> vuetify.VCheckbox:
        """
        Creates a Vuetify VCheckbox component with the following default components:
        - dense: set to 'true' to minimize space usage.

        :param label: Display label
        """
        InputComponents._build_component(
            vuetify.VCheckbox, label, v_model_name, **kwargs
        )

    @staticmethod
    def combobox(
        label: str,
        v_model_name: Optional[str] = None,
        items: Optional[list] = None,
        **kwargs,
    ) -> vuetify.VCombobox:
        """
        Creates a Vuetify VCombobox component with the following default components:
        - dense: set to 'true' to minimize space usage.

        :param label: Display label
        """
        InputComponents._build_component(
            vuetify.VCombobox, label, v_model_name, items=items, **kwargs
        )
