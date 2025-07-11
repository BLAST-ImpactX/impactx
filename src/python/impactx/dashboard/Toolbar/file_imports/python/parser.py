"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import state
from .helper import DashboardParserHelper

state.import_file = False
state.import_file_details = None
state.import_file_error = False
state.importing_file = False
state.imported_file_name = None


class DashboardParser:
    """
    Provides functionality to import ImpactX simulation files
    to the dashboard and auto-populate the UI with their configurations.
    """

    @staticmethod
    def reset_importing_states():
        """
        Resets import related states to default.
        """

        state.import_file_error = None
        state.import_file_details = None
        state.import_file = None
        state.importing_file = False
        state.imported_file_name = None

    @staticmethod
    def file_details(file) -> None:
        """
        Displays the size of the imported simulation file.

        :param file: ImpactX simulation file uploaded by the user.
        """

        file_size_in_bytes = file["size"]
        size_str = ""

        if file_size_in_bytes < 1024:
            size_str = f"{file_size_in_bytes} B"
        elif file_size_in_bytes < 1024 * 1024:
            size_str = f"{file_size_in_bytes / 1024:.1f} KB"

        state.imported_file_name = file["name"]
        state.import_file_details = f"({size_str}) {state.imported_file_name}"

    @staticmethod
    def parse_impactx_simulation_file(file) -> None:
        """
        Parses ImpactX simulation file contents.

        :param file: ImpactX simulation file uploaded by the user.
        """

        file_content = DashboardParserHelper.import_file_content(file, state)

        single_input_contents = DashboardParserHelper.parse_single_inputs(file_content)
        list_input_contents = DashboardParserHelper.parse_list_inputs(file_content)
        distribution_contents = DashboardParserHelper.parse_distribution(file_content)
        lattice_element_contents = DashboardParserHelper.parse_lattice_elements(
            file_content
        )
        used_vars = lattice_element_contents.pop("used_lattice_variables", set())

        variable_contents = DashboardParserHelper.parse_variables(
            file_content, used_vars
        )
        parsed_values_dictionary = {
            **single_input_contents,
            **list_input_contents,
            **distribution_contents,
            **lattice_element_contents,
            "variables": variable_contents,
        }

        return parsed_values_dictionary
