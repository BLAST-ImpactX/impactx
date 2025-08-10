"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from impactx import distribution

from ... import ctrl, state, vuetify
from ...Input.components import CardBase, CardComponents, InputComponents
from .. import DashboardDefaults
from ..defaults_helper import InputDefaultsHelper
from ..utils import GeneralFunctions
from ..validation import DashboardValidation, errors_tracker
from .utils import DistributionFunctions

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------

DISTRIBUTION_MODULE_NAME = distribution
DISTRIBUTION_PARAMETERS_AND_DEFAULTS = (
    InputDefaultsHelper.class_parameters_with_defaults(DISTRIBUTION_MODULE_NAME)
)

state.selected_distribution_parameters = {}
state.distribution_type_disable = False


def populate_distribution_parameters():
    """
    Called when `state.distribution_type` changes.
    Populates distribution parameters based on the current `state.distribution`.
    """
    params = {}
    is_twiss = state.distribution_type == "Twiss"

    # Gather necessary data
    if is_twiss:
        param_data = DistributionFunctions.get_twiss_data()
    else:
        # data for quadratic (impactX native)
        param_data = DISTRIBUTION_PARAMETERS_AND_DEFAULTS.get(state.distribution, [])

    # Populate the UI
    for param_name, default_value, default_type in param_data:
        error_message = DashboardValidation.validate(
            param_name, default_value, category="distribution"
        )
        units = DistributionFunctions.get_distribution_units(param_name)
        step = GeneralFunctions.get_default(param_name, "steps")

        params[param_name] = {
            "value": default_value,
            "type": default_type,
            "error_message": error_message,
            "units": units,
            "step": step,
        }

    state.selected_distribution_parameters = params
    errors_tracker.update_simulation_validation_status()
    return params


# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


@state.change("distribution")
def on_distribution_name_change(distribution, **kwargs):
    if distribution == "Thermal" or distribution == "Empty":
        state.distribution_type = ""
        state.distribution_type_disable = True
        state.dirty("distribution_type")
    else:
        type_list_default = DashboardDefaults.LISTS["distribution_type_list"]
        type_default = DashboardDefaults.DISTRIBUTION_PARAMETERS["distribution_type"]

        if state.distribution_type not in type_list_default:
            state.distribution_type = type_default

        state.distribution_type_disable = False


@state.change("distribution_type")
def on_distribution_type_change(**kwargs):
    populate_distribution_parameters()


@ctrl.add("update_distribution_parameter")
def on_distribution_parameter_change(name: str, input: str):
    numeric_input = GeneralFunctions.convert_to_numeric(input)
    error_message = DashboardValidation.validate(
        name, numeric_input, category="distribution"
    )

    parameter = state.selected_distribution_parameters.get(name)
    if parameter:
        parameter["value"] = numeric_input
        parameter["error_message"] = error_message
        errors_tracker.update_simulation_validation_status()
        state.dirty("selected_distribution_parameters")


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class DistributionConfiguration(CardBase):
    """
    User-Input section for beam distribution.
    """

    HEADER_NAME = "Distribution Parameters"

    def __init__(self):
        super().__init__()

    def card_content(self):
        """
        Creates UI content for beam distribution.
        """
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(self.HEADER_NAME)
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=6):
                        InputComponents.select(
                            label="Select Distribution",
                            v_model_name="distribution",
                        )
                    with vuetify.VCol(cols=6):
                        InputComponents.select(
                            label="Type",
                            v_model_name="distribution_type",
                            disabled=("distribution_type_disable",),
                        )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(
                        v_for="(parameter, parameter_name) in selected_distribution_parameters",
                        cols=4,
                    ):
                        with vuetify.VTooltip(
                            location="top",
                            text=("all_tooltips[parameter_name]",),
                        ):
                            with vuetify.Template(v_slot_activator="{ props }"):
                                vuetify.VTextField(
                                    label=("parameter_name",),
                                    v_model=("parameter.value",),
                                    id=("parameter_name",),
                                    suffix=("parameter.units",),
                                    update_modelValue=(
                                        ctrl.update_distribution_parameter,
                                        "[parameter_name, $event]",
                                    ),
                                    error_messages=("parameter.error_message",),
                                    type="number",
                                    step=("parameter.step",),
                                    __properties=["step"],
                                    density="compact",
                                    variant="underlined",
                                    hide_details="auto",
                                    v_bind="props",
                                )
