"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from impactx import elements

from .. import (
    CardComponents,
    InputComponents,
    NavigationComponents,
    generalFunctions,
    setup_server,
    vuetify,
)

server, state, ctrl = setup_server()

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
                "parameter_default_value": parameter[1],
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


def update_latticeElement_parameters(
    index, parameterName, parameterValue, parameterErrorMessage
):
    """
    Updates parameter value and includes error message if user input is not valid
    """

    for param in state.selected_lattice_list[index]["parameters"]:
        if param["parameter_name"] == parameterName:
            param["parameter_default_value"] = parameterValue
            param["parameter_error_message"] = parameterErrorMessage

    generalFunctions.update_simulation_validation_status()
    state.dirty("selected_lattice_list")


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
                    f"'{parameter['parameter_default_value']}'"
                )
            else:
                parameter_input[parameter["parameter_name"]] = parameter[
                    "parameter_default_value"
                ]
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


@ctrl.add("updateLatticeElementParameters")
def on_lattice_element_parameter_change(
    index, parameter_name, parameter_value, parameter_type
):
    parameter_value, input_type = generalFunctions.determine_input_type(parameter_value)
    error_message = generalFunctions.validate_against(parameter_value, parameter_type)

    update_latticeElement_parameters(
        index, parameter_name, parameter_value, error_message
    )


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


class LatticeConfiguration:
    """
    User-Input section for configuring lattice elements.
    """

    @staticmethod
    def card():
        """
        Creates UI content for lattice configuration.
        """

        with vuetify.VDialog(v_model=("showDialog", False)):
            LatticeConfiguration.dialog_configuration_list()

        with vuetify.VDialog(
            v_model=("lattice_configuration_dialog_settings", False), width="500px"
        ):
            LatticeConfiguration.dialog_settings()

        with vuetify.VCard(style="width: 696px;"):
            CardComponents.input_header("Lattice Configuration")
            with vuetify.VCardText():
                with vuetify.VRow(align="center", no_gutters=True):
                    with vuetify.VCol(cols=10):
                        vuetify.VCombobox(
                            label="Select Accelerator Lattice",
                            v_model=("selected_lattice", None),
                            items=("listOfLatticeElements",),
                            error_messages=("isSelectedLatticeListEmpty",),
                            dense=True,
                            classes="mr-2 pt-6",
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VBtn(
                            "ADD",
                            color="primary",
                            dense=True,
                            classes="mr-2",
                            click=ctrl.add_latticeElement,
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VIcon(
                            "mdi-cog",
                            click="lattice_configuration_dialog_settings = true",
                        )
                with vuetify.VRow():
                    with vuetify.VCol():
                        with vuetify.VCard(
                            style="height: 300px; width: 700px; overflow-y: auto;"
                        ):
                            with vuetify.VCardTitle(
                                "Elements", classes="text-subtitle-2 pa-3"
                            ):
                                vuetify.VSpacer()
                                vuetify.VIcon(
                                    "mdi-arrow-expand",
                                    color="primary",
                                    click="showDialog = true",
                                )
                            vuetify.VDivider()
                            LatticeConfiguration.configuration_list()

    # -----------------------------------------------------------------------------
    # Dialogs
    # -----------------------------------------------------------------------------

    @staticmethod
    def dialog_configuration_list():
        """
        Displays the configuration for lattice elements
        shown in a dialog box.
        """

        with vuetify.VCard():
            with vuetify.VCardTitle("Elements", classes="headline d-flex align-center"):
                vuetify.VSpacer()
                with vuetify.VBtn(icon=True, click="showDialog = false"):
                    vuetify.VIcon("mdi-close")
            vuetify.VDivider()
            LatticeConfiguration.configuration_list()

    @staticmethod
    def dialog_settings():
        """
        Provides controls for lattice element configuration,
        allowing dashboard users to define parameter defaults.
        """
        dialog_name = "lattice_configuration_dialog_tab_settings"

        NavigationComponents.create_dialog_tabs(dialog_name, 1, ["Defaults"])
        with vuetify.VTabsItems(v_model=(dialog_name, 0)):
            with vuetify.VTabItem():
                with vuetify.VCardText():
                    with vuetify.VRow():
                        with vuetify.VCol(cols=3):
                            InputComponents.text_field(
                                label="nslice",
                                v_model_name="nslice",
                                change=(
                                    ctrl.nsliceDefaultChange,
                                    "['nslice', $event]",
                                ),
                            )

    # -----------------------------------------------------------------------------
    # lattice_configuration_lsit
    # -----------------------------------------------------------------------------

    @staticmethod
    def configuration_list():
        """
        Displays the configuration for lattice elements.
        """
        with vuetify.VContainer(fluid=True):
            with vuetify.VRow(
                v_for="(latticeElement, index) in selected_lattice_list",
                align="center",
                no_gutters=True,
                style="min-width: 1500px;",
            ):
                with vuetify.VCol(cols="auto", classes="pa-2"):
                    vuetify.VIcon(
                        "mdi-menu-up",
                        click=(ctrl.move_latticeElementIndex_up, "[index]"),
                    )
                    vuetify.VIcon(
                        "mdi-menu-down",
                        click=(ctrl.move_latticeElementIndex_down, "[index]"),
                    )
                    vuetify.VIcon(
                        "mdi-delete",
                        click=(ctrl.deleteLatticeElement, "[index]"),
                    )
                    vuetify.VChip(
                        v_text=("latticeElement.name",),
                        dense=True,
                        classes="mr-2",
                        style="justify-content: center",
                    )
                with vuetify.VCol(
                    v_for="(parameter, parameterIndex) in latticeElement.parameters",
                    cols="auto",
                    classes="pa-2",
                ):
                    vuetify.VTextField(
                        label=("parameter.parameter_name",),
                        v_model=("parameter.parameter_default_value",),
                        change=(
                            ctrl.updateLatticeElementParameters,
                            "[index, parameter.parameter_name, $event, parameter.parameter_type]",
                        ),
                        error_messages=("parameter.parameter_error_message",),
                        dense=True,
                        style="width: 100px;",
                    )
