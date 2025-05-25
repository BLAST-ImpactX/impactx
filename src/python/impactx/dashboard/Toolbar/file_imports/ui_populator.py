"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server
from ...Input.distributionParameters.distributionMain import (
    on_distribution_parameter_change,
    populate_distribution_parameters,
)
from ...Input.latticeConfiguration.latticeMain import (
    add_lattice_element,
    on_lattice_element_parameter_change,
)
from .python.parser import DashboardParser
server, state, ctrl = setup_server()


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

    imported_distribution_data = imported_data["distribution"]["parameters"].items()
    imported_lattice_data = imported_data["lattice_elements"]
    non_state_inputs = ["distribution", "lattice_elements"]

    # Update state inputs (inputParameters, Space Charge, CSR, ISR)
    for input_name, input_value in imported_data.items():
        if hasattr(state, input_name) and input_name not in non_state_inputs:
            setattr(state, input_name, input_value)

    # Update distribution inputs
    if imported_distribution_data:
        state.distribution = imported_data["distribution"]["name"]
        state.distribution_type = imported_data["distribution"]["type"]
        state.flush()
        populate_distribution_parameters()

        for (
            distr_parameter_name,
            distr_parameter_value,
        ) in imported_distribution_data:
            on_distribution_parameter_change(
                distr_parameter_name, distr_parameter_value, "float"
            )

    # Update lattice elements
    state.selected_lattice_list = []

    for lattice_element_index, element in enumerate(imported_lattice_data):
        parsed_element = element["element"]
        parsed_parameters = element["parameters"]

        state.selected_lattice = parsed_element
        add_lattice_element()

        lattice_list_parameters = state.selected_lattice_list[
            lattice_element_index
        ]["parameters"]

        for (
            parsed_parameter_name,
            parsed_parameter_value,
        ) in parsed_parameters.items():
            parameter_type = None

            for parameter_info in lattice_list_parameters:
                parameter_info_name = parameter_info["parameter_name"]
                if parameter_info_name == parsed_parameter_name:
                    parameter_type = parameter_info["parameter_type"]
                    break

            if parameter_type:
                on_lattice_element_parameter_change(
                    lattice_element_index,
                    parsed_parameter_name,
                    parsed_parameter_value,
                    parameter_type,
                )