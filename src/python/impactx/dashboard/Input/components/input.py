from typing import Optional

from ... import vuetify
from ..shared import DROPDOWN_INPUTS
from ..utils import GeneralFunctions


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
            v_model_name = GeneralFunctions.normalize_for_v_model(label)

        if "items" in component_kwargs and component_kwargs["items"] is None:
            component_kwargs["items"] = (
                GeneralFunctions.get_default(f"{v_model_name}_list", "default_values"),
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
            location="top",
            text=(f"all_tooltips['{v_model_name}']",),
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

        if v_model_name is None:
            v_model_name = GeneralFunctions.normalize_for_v_model(label)

        DROPDOWN_INPUTS.add(v_model_name)

        InputComponents._build_component(
            vuetify.VSelect, label, v_model_name, items=items, **kwargs
        )

    @staticmethod
    def text_field(
        label: str, v_model_name: Optional[str] = None, input_type="number", **kwargs
    ) -> vuetify.VTextField:
        """
        Creates a Vuetify VTextField component with default properties.

        :param label: The label to display.
        :param v_model_name: Optional binding name for v_model. Defaults to a lowercase version
        :param input_type: The HTML input type (e.g., "text", "number", "email", etc.). Defaults to "number".
        of the label with spaces replaced by underscores if not provided.
        :param kwargs: Additional keyword arguments to pass to the component.
        """

        computed_v_model = (
            v_model_name
            if v_model_name is not None
            else GeneralFunctions.normalize_for_v_model(label)
        )
        InputComponents._build_component(
            vuetify.VTextField,
            label,
            v_model_name,
            error_messages=(f"{computed_v_model}_error_message", []),
            type=input_type,
            step=GeneralFunctions.get_default(computed_v_model, "steps"),
            suffix=GeneralFunctions.get_default(computed_v_model, "units"),
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
            vuetify.VSwitch, label, v_model_name, color="primary", **kwargs
        )

    @staticmethod
    def autocomplete(
        label: str,
        v_model_name: Optional[str] = None,
        items: Optional[list] = None,
        **kwargs,
    ) -> vuetify.VAutocomplete:
        """
        Creates a Vuetify VAutocomplete component with default properties.

        :param label: The label to display.
        :param v_model_name: Optional binding name for v_model. Defaults to a lowercase version
        of the label with spaces replaced by underscores if not provided.
        :param items: Optional list of items. If None, default items from defaults.py will be used.
        :param kwargs: Additional keyword arguments to pass to the component.
        """

        InputComponents._build_component(
            vuetify.VAutocomplete, label, v_model_name, items=items, **kwargs
        )
