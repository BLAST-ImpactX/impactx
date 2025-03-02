from typing import Optional

from ... import setup_server, vuetify
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
            "density": "compact",
            "variant": "underlined",
            "v_bind": "props",
            "hide_details": "auto",
        }
        props = {**common_props, **component_kwargs}

        with vuetify.VTooltip(
            location="bottom", text=TooltipDefaults.TOOLTIP.get(v_model_name)
        ):
            with vuetify.Template(v_slot_activator="{ props }"):
                vuetify_component(**props)

    @staticmethod
    def select(
        label: str,
        v_model_name: Optional[str] = None,
        items: Optional[list] = None,
        **kwargs,
    ) -> vuetify.VSelect:
        """
        Creates a Vuetify Vselect component with default properties.

        :param label: The label to display.
        :param v_model_name: Optional binding name for v_model. Defaults to a lowercase version
        of the label with spaces replaced by underscores if not provided.
        :param items: Optional list of items. If None, default items from defaults.py will be used.
        :param kwargs: Additional keyword arguments to pass to the component.
        """

        InputComponents._build_component(
            vuetify.VSelect, label, v_model_name, items=items, **kwargs
        )

    @staticmethod
    def text_field(
        label: str, v_model_name: Optional[str] = None, **kwargs
    ) -> vuetify.VTextField:
        """
        Creates a Vuetify VTextField component with default properties.

        :param label: The label to display.
        :param v_model_name: Optional binding name for v_model. Defaults to a lowercase version
        of the label with spaces replaced by underscores if not provided.
        :param kwargs: Additional keyword arguments to pass to the component.
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
        Creates a Vuetify VCheckbox component with default properties.

        :param label: The label to display.
        :param v_model_name: Optional binding name for v_model. Defaults to a lowercase version
        of the label with spaces replaced by underscores if not provided.
        :param kwargs: Additional keyword arguments to pass to the component.
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
        Creates a Vuetify VCombobox component with default properties.

        :param label: The label to display.
        :param v_model_name: Optional binding name for v_model. Defaults to a lowercase version
        of the label with spaces replaced by underscores if not provided.
        :param items: Optional list of items. If None, default items from defaults.py will be used.
        :param kwargs: Additional keyword arguments to pass to the component.
        """

        InputComponents._build_component(
            vuetify.VCombobox, label, v_model_name, items=items, **kwargs
        )
