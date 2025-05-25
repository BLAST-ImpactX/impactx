"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server
from .ui_populator import populate_impactx_simulation_file_to_ui
from .parser_helper import DashboardParserHelper

server, state, ctrl = setup_server()

state.imported_file_name = None


class ToolbarImport:
    @state.change("import_file")
    def on_import_file_change(import_file, **kwargs):
        if import_file:
            try:
                state.importing_file = True
                DashboardParser.file_details(import_file)
                parsed_data = DashboardParser.parse_impactx_simulation_file(import_file)
                populate_impactx_simulation_file_to_ui(parsed_data)
            except Exception:
                state.import_file_error = True
                state.import_file_error_message = "Unable to parse"
            finally:
                state.importing_file = False

    @staticmethod
    def reset_importing_states():
        """
        Resets import related states to default.
        """

        state.import_file_error = None
        state.import_file_details = None
        state.import_file = None
        state.importing_file = False


class DashboardParser:
    """
    Provides functionality to import ImpactX simulation files
    to the dashboard and auto-populate the UI with their configurations.
    """

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

        parsed_values_dictionary = {
            **single_input_contents,
            **list_input_contents,
            **distribution_contents,
            **lattice_element_contents,
        }

        return parsed_values_dictionary