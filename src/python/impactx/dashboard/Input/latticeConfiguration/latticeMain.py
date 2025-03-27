"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from impactx import elements

from ... import setup_server, vuetify
from .. import (
    CardBase,
    CardComponents,
    InputComponents,
    generalFunctions,
)
from . import LatticeConfigurationHelper, LatticeVariableHandler

server, state, ctrl = setup_server()
state.lattice_params_bound_or_pending_variable = {}

# -----------------------------------------------------------------------------
# Helpful
# -----------------------------------------------------------------------------

LATTICE_ELEMENTS_MODULE_NAME = elements

state.listOfLatticeElements = generalFunctions.select_classes(
    LATTICE_ELEMENTS_MODULE_NAME
)
state.listOfLatticeElementParametersAndDefault = (
    generalFunctions.class_parameters_with_defaults(LATTICE_ELEMENTS_MODULE_NAME)
)

# -----------------------------------------------------------------------------
# Default
# -----------------------------------------------------------------------------

state.selected_lattice_list = []
state.nslice = ""

# -----------------------------------------------------------------------------
# Main Functions
# -----------------------------------------------------------------------------


def add_lattice_element():
    """
    Adds the selected lattice element to the list of selected
    lattice elements along with its default parameters.
    :return: dictionary representing the added lattice element with its parameters.
    """

    selected_lattice = state.selected_lattice
    selected_lattice_parameters = state.listOfLatticeElementParametersAndDefault.get(
        selected_lattice, []
    )

    selected_lattice_element = {
        "name": selected_lattice,
        "parameters": [
            {
                "parameter_name": parameter[0],
                "ui_value": parameter[1],
                "sim_value": parameter[1],
                "parameter_type": parameter[2],
                "parameter_error_message": generalFunctions.validate_against(
                    parameter[1], parameter[2]
                ),
            }
            for parameter in selected_lattice_parameters
        ],
    }

    state.selected_lattice_list.append(selected_lattice_element)
    generalFunctions.update_simulation_validation_status()
    return selected_lattice_element


# -----------------------------------------------------------------------------
# Write to file functions
# -----------------------------------------------------------------------------


def parameter_input_checker_for_lattice(latticeElement):
    """
    Helper function to check if user input is valid.
    :return: A dictionary with parameter names as keys and their validated values.
    """

    parameter_input = {}
    for parameter in latticeElement["parameters"]:
        if parameter["parameter_error_message"] == []:
            if parameter["parameter_type"] == "str":
                parameter_input[parameter["parameter_name"]] = (
                    f"'{parameter['sim_value']}'"
                )
            else:
                parameter_input[parameter["parameter_name"]] = parameter["sim_value"]
        else:
            parameter_input[parameter["parameter_name"]] = 0

    return parameter_input


def lattice_elements():
    """
    Writes user input for lattice element parameters parameters in suitable format for simulation code.
    :return: A list in the suitable format.
    """

    elements_list = []
    for latticeElement in state.selected_lattice_list:
        latticeElement_name = latticeElement["name"]
        parameters = parameter_input_checker_for_lattice(latticeElement)

        param_values = ", ".join(f"{value}" for value in parameters.values())
        elements_list.append(eval(f"elements.{latticeElement_name}({param_values})"))

    return elements_list


# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


@state.change("selected_lattice_list")
def on_selected_lattice_list_change(selected_lattice_list, **kwargs):
    if selected_lattice_list == []:
        state.isSelectedLatticeListEmpty = "Please select a lattice element"
        generalFunctions.update_simulation_validation_status()
    else:
        state.isSelectedLatticeListEmpty = ""


@state.change("selected_lattice")
def on_lattice_element_name_change(selected_lattice, **kwargs):
    return


@ctrl.add("add_latticeElement")
def on_add_lattice_element_click():
    selected_lattice = state.selected_lattice

    if selected_lattice not in state.listOfLatticeElements:
        state.isSelectedLatticeListEmpty = (
            f"Lattice element '{selected_lattice}' does not exist."
        )
    else:
        add_lattice_element()
        state.dirty("selected_lattice_list")


def process_if_variable(index, parameter_name, ui_value, parameter_type):
    """
    If the updated lattice parameter value uses or potentially uses a variable, this
    function returns the simulation value by lookup, adds the element to a dictionary
    which contains current or potential variables, and also returns true or false if it
    is a current or potential variable.

    :param index: The index of the lattice element in the lattice list config.
    :param parameter_name: The specific lattice element parameter name.
    :param ui_value: The value present on the UI end..
    :param parameter_type: The lattice element parameters type.
    """

    lattice_variable, variable_index = LatticeVariableHandler.determine_if_existing_variable(
        ui_value
    )
    potentially_lattice_variable = (
        LatticeVariableHandler.element_potentially_using_element(ui_value, index)
    )

    if lattice_variable or potentially_lattice_variable:
        if lattice_variable and variable_index is not None:
            sim_value = state.variables[variable_index]["value"]
        else:
            sim_value = ui_value

        binding = {
            "index": index,
            "parameter_name": parameter_name,
            "ui_value": ui_value,
            "parameter_type": parameter_type,
            "variable_index": variable_index,
        }
    else:
        sim_value, _ = generalFunctions.determine_input_type(ui_value)
        binding = None

    return sim_value, binding


@ctrl.add("updateLatticeElementParameters")
def on_lattice_element_parameter_change(
    index, parameter_name, ui_value, parameter_type
):
    sim_value, bounded_or_pending_variable = process_if_variable(
        index, parameter_name, ui_value, parameter_type
    )

    if bounded_or_pending_variable is not None:
        state.lattice_params_bound_or_pending_variable[(index, parameter_name)] = (
            bounded_or_pending_variable
        )
    else:
        state.lattice_params_bound_or_pending_variable.pop(
            (index, parameter_name), None
        )

    error_message = generalFunctions.validate_against(sim_value, parameter_type)

    for param in state.selected_lattice_list[index]["parameters"]:
        if param["parameter_name"] == parameter_name:
            param["ui_value"] = ui_value
            param["sim_value"] = sim_value
            param["parameter_error_message"] = error_message

    print(
        f"Updated lattice parameter '{parameter_name}' at index {index}: ui_value={ui_value}, sim_value={sim_value}, error={error_message}"
    )
    print(
        f"# of elements using/potentially using variables: {len(state.lattice_params_bound_or_pending_variable)}"
    )

    generalFunctions.update_simulation_validation_status()
    state.dirty("selected_lattice_list")


@ctrl.add("deleteLatticeElement")
def on_delete_LatticeElement_click(index):
    state.selected_lattice_list.pop(index)
    state.dirty("selected_lattice_list")


@ctrl.add("move_latticeElementIndex_up")
def on_move_latticeElementIndex_up_click(index):
    if index > 0:
        state.selected_lattice_list[index], state.selected_lattice_list[index - 1] = (
            state.selected_lattice_list[index - 1],
            state.selected_lattice_list[index],
        )
        state.dirty("selected_lattice_list")


@ctrl.add("move_latticeElementIndex_down")
def on_move_latticeElementIndex_down_click(index):
    if index < len(state.selected_lattice_list) - 1:
        state.selected_lattice_list[index], state.selected_lattice_list[index + 1] = (
            state.selected_lattice_list[index + 1],
            state.selected_lattice_list[index],
        )
        state.dirty("selected_lattice_list")


@ctrl.add("nsliceDefaultChange")
def update_default_value(parameter_name, new_value):
    data = generalFunctions.class_parameters_with_defaults(elements)

    for key, parameters in data.items():
        for i, param in enumerate(parameters):
            if param[0] == parameter_name:
                parameters[i] = (param[0], new_value, param[2])

    state.listOfLatticeElementParametersAndDefault = data


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class LatticeConfiguration(CardBase):
    HEADER_NAME = "Lattice Configuration"

    def __init__(self):
        super().__init__()

    def init_settings_dialog(self):
        with vuetify.VDialog(
            v_model=("lattice_configuration_dialog_settings", False), width="500px"
        ):
            LatticeVariableHandler.dialog_settings()

    def card_content(self):
        self.init_settings_dialog()
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(
                self.HEADER_NAME,
                additional_components={
                    "end": LatticeConfigurationHelper.settings,
                },
            )
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=True):
                        InputComponents.combobox(
                            label="Select Accelerator Lattice",
                            v_model_name="selected_lattice",
                            items=("listOfLatticeElements",),
                            error_messages=("isSelectedLatticeListEmpty",),
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VBtn(
                            "ADD",
                            color="primary",
                            dense=True,
                            click=ctrl.add_latticeElement,
                        )
                with vuetify.VRow(
                    **self.ROW_STYLE,
                    v_for="(latticeElement, index) in selected_lattice_list",
                    align="center",
                    style="flex-wrap: nowrap;",
                ):
                    with vuetify.VCol(cols="auto"):
                        LatticeConfigurationHelper.move_element_up()
                        LatticeConfigurationHelper.move_element_down()
                        LatticeConfigurationHelper.delete_element()
                    with vuetify.VCol(cols="auto"):
                        vuetify.VChip(
                            text=("latticeElement.name",),
                            style="justify-content: center",
                        )
                    with vuetify.VCol(
                        v_for="(parameter, parameterIndex) in latticeElement.parameters",
                        cols="auto",
                    ):
                        vuetify.VTextField(
                            label=("parameter.parameter_name",),
                            v_model=("parameter.ui_value",),
                            update_modelValue=(
                                ctrl.updateLatticeElementParameters,
                                "[index, parameter.parameter_name, $event, parameter.parameter_type]",
                            ),
                            error_messages=("parameter.parameter_error_message",),
                            density="comfortable",
                            variant="underlined",
                            style="width: 100px;",
                        )
