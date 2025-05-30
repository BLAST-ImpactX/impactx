"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from impactx import distribution

from ... import setup_server, vuetify
from .. import (
    CardBase,
    CardComponents,
    DashboardDefaults,
    InputComponents,
    generalFunctions,
)
from . import DistributionFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------

DISTRIBUTION_MODULE_NAME = distribution
DISTRIBUTION_LIST = generalFunctions.select_classes(DISTRIBUTION_MODULE_NAME)
DISTRIBUTION_PARAMETERS_AND_DEFAULTS = generalFunctions.class_parameters_with_defaults(
    DISTRIBUTION_MODULE_NAME
)

state.selected_distribution_parameters = {}
state.distribution_type_disable = False


def populate_distribution_parameters():
    """
    Populates distribution parameters based on the selected distribution.
    """
    params = {}
    param_data = []
    is_twiss = state.distribution_type == "Twiss"

    # Gather necessary data
    if is_twiss:
        param_data = DistributionFunctions.get_twiss_data()
    else:
        # data for quadratic (impactX native)
        param_data = DISTRIBUTION_PARAMETERS_AND_DEFAULTS.get(state.distribution, [])

    # Populate the UI
    for param_name, param_value, param_type in param_data:
        error_message = generalFunctions.validate_against(param_value, param_type)
        units = DistributionFunctions.get_distribution_units(param_name)
        step = generalFunctions.get_default(param_name, "steps")

        params[param_name] = {
            "parameter_default_value": param_value,
            "parameter_type": param_type,
            "parameter_error_message": error_message,
            "parameter_units": units,
            "parameter_step": step,
        }

    state.selected_distribution_parameters = params
    generalFunctions.update_simulation_validation_status()
    return params


# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


@state.change("distribution")
def on_distribution_name_change(distribution, **kwargs):
    if state.importing_file:
        return

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
    if state.importing_file:
        return
    populate_distribution_parameters()


@ctrl.add("update_distribution_parameter")
def on_distribution_parameter_change(parameter_name, parameter_value, parameter_type):
    parameter_value = generalFunctions.convert_to_numeric(parameter_value)
    lookup_name = "lambda" if "lambda" in parameter_name else parameter_name
    conditions = generalFunctions.get_default(lookup_name, "validation_condition")
    error_message = generalFunctions.validate_against(
        parameter_value, parameter_type, additional_conditions=conditions
    )

    if parameter_name in state.selected_distribution_parameters:
        state.selected_distribution_parameters[parameter_name][
            "parameter_default_value"
        ] = parameter_value
        state.selected_distribution_parameters[parameter_name][
            "parameter_error_message"
        ] = error_message

    generalFunctions.update_simulation_validation_status()
    state.dirty("selected_distribution_parameters")


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class DistributionParameters(CardBase):
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
                            items=(DISTRIBUTION_LIST,),
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
                                    v_model=("parameter.parameter_default_value",),
                                    suffix=("parameter.parameter_units",),
                                    update_modelValue=(
                                        ctrl.update_distribution_parameter,
                                        "[parameter_name, $event, parameter.parameter_type]",
                                    ),
                                    error_messages=(
                                        "parameter.parameter_error_message",
                                    ),
                                    type="number",
                                    step=("parameter.parameter_step",),
                                    __properties=["step"],
                                    density="compact",
                                    variant="underlined",
                                    hide_details="auto",
                                    v_bind="props",
                                )
