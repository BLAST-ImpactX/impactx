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

    @staticmethod
    def checkbox(
        label: str, v_model_name: Optional[str] = None, **kwargs
    ) -> vuetify.VCheckbox:
        """
        Creates a Vuetify VCheckbox component with the following default components:
        - dense: set to 'true' to minimize space usage.

        :param label: Display label
        """
        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        with vuetify.VTooltip(**TooltipDefaults.TOOLTIP_STYLE):
            with vuetify.Template(v_slot_activator="{ on, attrs }"):
                vuetify.VCheckbox(
                    label=label,
                    v_model=(v_model_name,),
                    dense=True,
                    v_on="on",
                    v_bind="attrs",
                    **kwargs,
                )
            html.Span(TooltipDefaults.TOOLTIP.get(v_model_name))

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
        if v_model_name is None:
            v_model_name = label.lower().replace(" ", "_")

        with vuetify.VTooltip(**TooltipDefaults.TOOLTIP_STYLE):
            with vuetify.Template(v_slot_activator="{ on, attrs }"):
                vuetify.VCombobox(
                    label=label,
                    v_model=(v_model_name,),
                    items=items,
                    dense=True,
                    **kwargs,
                    v_on="on",
                    v_bind="attrs",
                )
            html.Span(TooltipDefaults.TOOLTIP.get(v_model_name))