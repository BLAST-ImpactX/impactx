"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import ctrl, state
from ...Input.lattice.ui import add_lattice_element
from ...Input.lattice.variable_handler import LatticeVariableHandler
from .python.parser import DashboardParser


@state.change("import_file")
def on_import_file_change(import_file, **kwargs):
    if import_file:
        try:
            state.importing_file = True
            DashboardParser.file_details(import_file)
            populate_impactx_simulation_file_to_ui(import_file)
        except Exception:
            state.import_file_error = True
            state.import_file_error_message = "Unable to parse"
        finally:
            state.importing_file = False


def populate_impactx_simulation_file_to_ui(file) -> None:
    """
    Auto fills the dashboard with parsed inputs.

    :param file: ImpactX simulation file uploaded by the user.
    """

    imported_data = DashboardParser.parse_impactx_simulation_file(file)

    imported_distribution_data = imported_data["distribution"]
    imported_lattice_data = imported_data["lattice_elements"]
    parsed_variables = imported_data["variables"]
    non_state_inputs = ["distribution", "lattice_elements", "variables"]

    # Update state inputs (simulation parameters, Space Charge, CSR, ISR)
    for input_name, input_value in imported_data.items():
        if hasattr(state, input_name) and input_name not in non_state_inputs:
            setattr(state, input_name, input_value)

    _populate_distribution_inputs_to_ui(imported_distribution_data)
    _populate_lattice_config_to_ui(imported_lattice_data)
    _populate_lattice_config_variables_to_ui(parsed_variables)


@staticmethod
def _populate_distribution_inputs_to_ui(parsed_data):
    # the below two calls do not call state.change("distribution") or state.change("distribution_type")
    # since they are both part of a nested state (ie. distribution_type=["Twiss","Quadratic"]).
    state.distribution = parsed_data["name"]
    state.distribution_type = parsed_data["type"]
    state.flush()  # force calls state.change("distribution") and state.change("distribution_type")

    for distr_param_name, distr_param_value in parsed_data["parameters"].items():
        ctrl.update_distribution_parameter(distr_param_name, distr_param_value)


@staticmethod
def _populate_lattice_config_to_ui(parsed_data):
    # Update lattice elements
    state.selected_lattice_list = []

    for lattice_element_index, element in enumerate(parsed_data):
        parsed_element = element["element"]
        parsed_parameters = element["parameters"]

        state.selected_lattice = parsed_element
        add_lattice_element()

        lattice_list_parameters = state.selected_lattice_list[lattice_element_index][
            "parameters"
        ]

        for parsed_param_name, parsed_param_value in parsed_parameters.items():
            parameter_type = None

            for parameter_info in lattice_list_parameters:
                parameter_info_name = parameter_info["parameter_name"]
                if parameter_info_name == parsed_param_name:
                    parameter_type = parameter_info["parameter_type"]
                    break

            if parameter_type:
                ctrl.updateLatticeElementParameters(
                    lattice_element_index,
                    parsed_param_name,
                    parsed_param_value,
                    parameter_type,
                )


@staticmethod
def _populate_lattice_config_variables_to_ui(parsed_data):
    # Remove default empty entry if it exists
    state.variables = [var for var in state.variables if var["name"]]

    for name, value in parsed_data.items():
        # Check if a variable with the same name already exists
        if not any(var["name"] == name for var in state.variables):
            state.variables.append({"name": name, "value": value, "error_message": ""})
    state.dirty("variables")
    LatticeVariableHandler.update_delete_availability()
