"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import inspect

from distribution_input_helpers import twiss
from trame.widgets import vuetify

from impactx import distribution

from ...Input.trameFunctions import TrameFunctions
from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions
from .distributionFunctions import DistributionFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------

DISTRIBUTIONS_MODULE_NAME = distribution

state.listOfDistributions = generalFunctions.select_classes(DISTRIBUTIONS_MODULE_NAME)
state.listOfDistributionsAndParametersAndDefault = (
    generalFunctions.class_parameters_with_defaults(DISTRIBUTIONS_MODULE_NAME)
)

# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------

state.selected_distribution = generalFunctions.get_default(
    "selected_distribution", "default_values"
)
state.selected_distribution_type = generalFunctions.get_default(
    "selected_distribution_type", "default_values"
)
state.selected_distribution_parameters = []
state.distributionTypeDisabled = False

# -----------------------------------------------------------------------------
# Main Functions
# -----------------------------------------------------------------------------


def populate_distribution_parameters(selected_distribution):
    """
    Populates distribution parameters based on the selected distribution.
    :param selected_distribution (str): The name of the selected distribution
    whose parameters need to be populated.
    """

    if state.selected_distribution_type == "Twiss":
        sig = inspect.signature(twiss)
        state.selected_distribution_parameters = [
            {
                "parameter_name": param.name,
                "parameter_default_value": param.default
                if param.default != param.empty
                else None,
                "parameter_type": "float",  # Hardcoding Twiss to 'float' type.
                "parameter_error_message": generalFunctions.validate_against(
                    param.default if param.default != param.empty else None, "float"
                ),
                "parameter_units": generalFunctions.get_default(param.name, "units")
                if "beta" in param.name or "emitt" in param.name
                else "",
                "parameter_step": generalFunctions.get_default(param.name, "steps"),
            }
            for param in sig.parameters.values()
        ]

    else:  # when type == 'Quadratic Form'
        selected_distribution_parameters = (
            state.listOfDistributionsAndParametersAndDefault.get(
                selected_distribution, []
            )
        )

        state.selected_distribution_parameters = [
            {
                "parameter_name": parameter[0],
                "parameter_default_value": parameter[1],
                "parameter_type": parameter[2],
                "parameter_error_message": generalFunctions.validate_against(
                    parameter[1], parameter[2]
                ),
                "parameter_units": "m"
                if "beta" in parameter[0] or "emitt" in parameter[0]
                else "",
                "parameter_step": generalFunctions.get_default(parameter[0], "steps"),
            }
            for parameter in selected_distribution_parameters
        ]

    generalFunctions.update_simulation_validation_status()
    return state.selected_distribution_parameters


def update_distribution_parameters(
    parameterName, parameterValue, parameterErrorMessage
):
    """
    Updates the value of a distribution parameter and its error message.

    :param parameterName (str): The name of the parameter to update.
    :param parameterValue: The new value for the parameter.
    :param parameterErrorMessage: The error message related to the parameter's value.
    """

    for param in state.selected_distribution_parameters:
        if param["parameter_name"] == parameterName:
            param["parameter_default_value"] = parameterValue
            param["parameter_error_message"] = parameterErrorMessage

    generalFunctions.update_simulation_validation_status()
    state.dirty("selected_distribution_parameters")


# -----------------------------------------------------------------------------
# Write to file functions
# -----------------------------------------------------------------------------


def distribution_parameters():
    """
    :return: An instance of the selected distribution class,
    initialized with the appropriate parameters provided by the user.
    """

    distribution_name = state.selected_distribution
    parameters = DistributionFunctions.convert_distribution_parameters_to_valid_type()

    if state.selected_distribution_type == "Twiss":
        twiss_params = twiss(**parameters)
        distr = getattr(distribution, distribution_name)(**twiss_params)
    else:
        distr = getattr(distribution, distribution_name)(**parameters)

    return distr


# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


@state.change("selected_distribution")
def on_distribution_name_change(selected_distribution, **kwargs):
    if selected_distribution == "Thermal":
        state.selected_distribution_type = "Quadratic Form"
        state.distributionTypeDisabled = True
        state.dirty("selected_distribution_type")
    else:
        state.distributionTypeDisabled = False
    populate_distribution_parameters(selected_distribution)


@state.change("selected_distribution_type")
def on_distribution_type_change(**kwargs):
    populate_distribution_parameters(state.selected_distribution)


@ctrl.add("updateDistributionParameters")
def on_distribution_parameter_change(parameter_name, parameter_value, parameter_type):
    parameter_value, input_type = generalFunctions.determine_input_type(parameter_value)
    error_message = generalFunctions.validate_against(parameter_value, parameter_type)

    update_distribution_parameters(parameter_name, parameter_value, error_message)


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class DistributionParameters:
    """
    User-Input section for beam distribution.
    """

    @staticmethod
    def card():
        """
        Creates UI content for beam distribution.
        """

        with vuetify.VCard(style="width: 340px; height: 300px"):
            with vuetify.VCardTitle("Distribution Parameters"):
                vuetify.VSpacer()
                TrameFunctions.create_refresh_button(
                    lambda: generalFunctions.reset_inputs("distribution")
                )
                vuetify.VIcon(
                    "mdi-information",
                    style="color: #00313C;",
                    click=lambda: generalFunctions.documentation("BeamDistributions"),
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow():
                    with vuetify.VCol(cols=6):
                        vuetify.VCombobox(
                            label="Select Distribution",
                            v_model=("selected_distribution",),
                            items=("listOfDistributions",),
                            dense=True,
                        )
                    with vuetify.VCol(cols=6):
                        vuetify.VSelect(
                            v_model=("selected_distribution_type",),
                            label="Type",
                            items=(
                                generalFunctions.get_default(
                                    "distribution_type_list", "default_values"
                                ),
                            ),
                            dense=True,
                            disabled=("distributionTypeDisabled",),
                        )
                with vuetify.VRow(classes="my-2"):
                    for i in range(3):
                        with vuetify.VCol(cols=4, classes="py-0"):
                            with vuetify.VRow(
                                v_for="(parameter, index) in selected_distribution_parameters"
                            ):
                                with vuetify.VCol(
                                    v_if=f"index % 3 == {i}", classes="py-1"
                                ):
                                    vuetify.VTextField(
                                        label=("parameter.parameter_name",),
                                        v_model=("parameter.parameter_default_value",),
                                        suffix=("parameter.parameter_units",),
                                        change=(
                                            ctrl.updateDistributionParameters,
                                            "[parameter.parameter_name, $event, parameter.parameter_type]",
                                        ),
                                        error_messages=(
                                            "parameter.parameter_error_message",
                                        ),
                                        type="number",
                                        step=("parameter.parameter_step",),
                                        __properties=["step"],
                                        dense=True,
                                    )
