"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from typing import Union

from .. import setup_server
from ..Toolbar.file_imports.python.parser import DashboardParser
from .defaults import DashboardDefaults

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------


class generalFunctions:
    @staticmethod
    def open_documentation(section_name):
        """
        Retrieves the documentation link with the provided section_name
        and opens the documentation sidebar on the dashoard.

        :param section_name: The name for the input section.
        """

        new_url = DashboardDefaults.DOCUMENTATION.get(section_name)
        if state.documentation_drawer_open and state.documentation_url == new_url:
            state.documentation_drawer_open = False
        else:
            state.documentation_url = new_url
            state.documentation_drawer_open = True

    @staticmethod
    def get_default(parameter, type):
        parameter_type_dictionary = getattr(DashboardDefaults, f"{type.upper()}", None)
        parameter_default = parameter_type_dictionary.get(parameter)

        if parameter_default is not None:
            return parameter_default

        parameter_name_base = parameter.partition("_")[0]
        return parameter_type_dictionary.get(parameter_name_base)

    # -----------------------------------------------------------------------------
    # Validation functions
    # -----------------------------------------------------------------------------

    @staticmethod
    def convert_to_numeric(input: str) -> Union[int, float]:
        """
        Converts string inputs to their appropriate numeric type.
        This method is needed since text fields inputs on the dashboard
        are inherently strings.

        It first tries to convert the value to int, then to float.
        If the conversion fails, returns None.
        Note that the function runs on every keystroke.
        For non-trivial input, e.g., '1e-2', the conversion
        fails silently until the full number is typed.

        :param input: The input to convert to a numeric type.
        :return: The input converted to a numeric type.
        """

        try:
            return int(input)
        except ValueError:
            try:
                return float(input)
            except ValueError:
                return None

    @staticmethod
    def reset_inputs(input_section):
        """
        Resets dashboard inputs to default values.

        :param input_section: The input section to reset.
        """

        possible_section_names = []
        for name in vars(DashboardDefaults):
            if name != "DEFAULT_VALUES" and name.isupper():
                possible_section_names.append(name)

        if input_section.upper() in possible_section_names:
            state.update(getattr(DashboardDefaults, input_section.upper()))

            if input_section == "distribution_parameters":
                state.dirty("distribution_type")
            elif input_section == "lattice_configuration":
                state.selected_lattice_list = []
                state.variables = [{"name": "", "value": "", "error_message": ""}]
                state.dirty("variables")
            elif input_section == "space_charge":
                state.dirty("max_level")

        elif input_section == "all":
            DashboardParser.reset_importing_states()
            state.update(DashboardDefaults.DEFAULT_VALUES)
            state.dirty("distribution_type")
            state.selected_lattice_list = []
            state.dirty("max_level")
            state.variables = [{"name": "", "value": "", "error_message": ""}]
            state.dirty("variables")
