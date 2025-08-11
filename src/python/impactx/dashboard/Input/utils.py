"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from typing import Union

from .. import state
from ..Toolbar.file_imports.python.parser import DashboardParser
from .defaults import DashboardDefaults


class GeneralFunctions:
    @staticmethod
    def normalize_for_v_model(name: str) -> str:
        """
        Normalizes a name for use as a v-model variable name.
        Converts to lowercase with spaces replaced by underscores.

        :param name: The name to normalize
        :return: Normalized v-model name
        """
        return name.lower().replace(" ", "_")

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
    def get_default(parameter: str, type: str) -> str | None:
        """
        Get the default value for a parameter by exact or base name match.

        Attempts full match first, then falls back to removing the last underscore suffix.

        :param parameter: Full parameter name (e.g., 'beta_x', 'blocking_factor_x').
        :param type: Parameter group name (e.g., 'simulation_parameters', 'csr').
        :return: Default value if found, else None.
        """
        parameter_type_dictionary = getattr(DashboardDefaults, f"{type.upper()}", None)
        if not parameter_type_dictionary:
            return None

        if parameter in parameter_type_dictionary:
            return parameter_type_dictionary[parameter]

        if "_" in parameter:
            parameter_name_base = "_".join(parameter.split("_")[:-1])
            return parameter_type_dictionary.get(parameter_name_base)

        return None

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

        if isinstance(input, (int, float)):
            return input

        try:
            return int(input)
        except (ValueError, TypeError):
            try:
                return float(input)
            except (ValueError, TypeError):
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

    @staticmethod
    def set_state_to_numeric(state_name: str) -> None:
        """
        Converts the value of a state variable to a numeric type (int or float)
        and updates the state in-place.

        :param state_name: The name of the state variable to convert and update.
        """
        current_input = getattr(state, state_name)
        numeric_input = GeneralFunctions.convert_to_numeric(current_input)
        setattr(state, state_name, numeric_input)
