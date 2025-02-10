"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import inspect
import re

from .. import setup_server
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
    def determine_input_type(value):
        """
        Determines the type of the input value.
        :param value: The input value whose type needs to be determined.
        :return: A tuple containing the value converted to its determined type and the type itself.
        """

        try:
            return int(value), int
        except ValueError:
            try:
                return float(value), float
            except ValueError:
                return value, str

    @staticmethod
    def validate_against(input_value, value_type, additional_conditions=None):
        """
        Validates the input value against the desired type and additional conditions.
        :param input_value: The value to validate.
        :param value_type: The desired type ('int', 'float', 'str').
        :param additional_conditions: A list of additional conditions to validate.
        :return: A list of error messages. An empty list if there are no errors.
        """
        errors = []
        value = None

        if input_value == "None":
            return errors

        # value_type checking
        if value_type == "int":
            if input_value is None:
                errors.append("Must be an integer")
            else:
                try:
                    value = int(input_value)
                except ValueError:
                    errors.append("Must be an integer")
        elif value_type == "float":
            if input_value is None:
                errors.append("Must be a float")
            else:
                try:
                    value = float(input_value)
                except ValueError:
                    errors.append("Must be a float")
        elif value_type == "str":
            if input_value is None:
                errors.append("Must be a string")
            else:
                value = str(input_value)
        else:
            errors.append("Unknown type")

        # addition_conditions checking
        if errors == [] and additional_conditions:
            for condition in additional_conditions:
                if condition == "non_zero" and value == 0:
                    errors.append("Must be non-zero.")
                if condition == "positive" and value < 0:
                    errors.append("Must be positive.")
                if condition == "negative" and value > 0:
                    errors.append("Must be negative.")

        return errors

    @staticmethod
    def update_simulation_validation_status():
        """
        Checks if any input fields are not provided with the correct input type.
        Updates the state to enable or disable the run simulation button.
        """

        error_details = []

        # Check for errors in distribution parameters
        for param in state.selected_distribution_parameters:
            if param["parameter_error_message"]:
                error_details.append(
                    f"{param['parameter_name']}: {param['parameter_error_message']}"
                )

        # Check for errors in lattice parameters
        for lattice in state.selected_lattice_list:
            for param in lattice["parameters"]:
                if param["parameter_error_message"]:
                    error_details.append(
                        f"Lattice {lattice['name']} - {param['parameter_name']}: {param['parameter_error_message']}"
                    )

        # Check for errors in input card
        if state.npart_error_message:
            error_details.append(f"Number of Particles: {state.npart_error_message}")
        if state.kin_energy_error_message:
            error_details.append(f"Kinetic Energy: {state.kin_energy_error_message}")
        if state.bunch_charge_C_error_message:
            error_details.append(f"Bunch Charge: {state.bunch_charge_C_error_message}")
        if state.charge_qe_error_message:
            error_details.append(
                f"Ref. Particle Charge: {state.charge_qe_error_message}"
            )
        if state.mass_MeV_error_message:
            error_details.append(f"Ref. Particle Mass: {state.mass_MeV_error_message}")

        if state.selected_lattice_list == []:
            error_details.append("LatticeListIsEmpty")

        # Check for errors in CSR parameters
        if state.csr_bins_error_message:
            error_details.append(f"CSR Bins: {state.csr_bins_error_message}")

        # Check for errors in Space Charge parameters
        if state.space_charge:
            # n_cell parameters
            for direction in ["x", "y", "z"]:
                n_cell_error = getattr(state, f"error_message_n_cell_{direction}")
                if n_cell_error:
                    error_details.append(f"n_cell_{direction}: {n_cell_error}")

            # Blocking factor parameters
            for direction in ["x", "y", "z"]:
                blocking_factor_error = getattr(
                    state, f"error_message_blocking_factor_{direction}"
                )
                if blocking_factor_error:
                    error_details.append(
                        f"blocking_factor_{direction}: {blocking_factor_error}"
                    )

            # Prob Relative Fields
            for index, field in enumerate(state.prob_relative_fields):
                if field["error_message"]:
                    error_details.append(
                        f"prob_relative[{index}]: {field['error_message']}"
                    )

        state.disableRunSimulationButton = bool(error_details)

    # -----------------------------------------------------------------------------
    # Class, parameter, default value, and default type retrievals
    # -----------------------------------------------------------------------------

    @staticmethod
    def find_classes(module_name):
        """
        Returns a list of all classes in the given module.
        :param module_name: The module to inspect.
        :return: A list of tuples containing class names.
        """

        results = []
        for name in dir(module_name):
            attr = getattr(module_name, name)
            if inspect.isclass(attr):
                results.append((name, attr))
        return results

    @staticmethod
    def find_init_docstring_for_classes(classes):
        """
        Retrieves the __init__ docstring of the given classes.
        :param classes: A list of typles containing class names.
        :return: A dictionary with class names as keys and their __init__ docstrings as values.
        """

        if not isinstance(classes, (list, tuple)):
            raise TypeError("The 'classes' argument must be a list or tuple.")

        docstrings = {}
        for name, cls in classes:
            init_method = getattr(cls, "__init__", None)
            if init_method:
                docstring = cls.__init__.__doc__
                docstrings[name] = docstring
        return docstrings

    @staticmethod
    def extract_parameters(docstring):
        """
        Parses specific information from docstrings.
        Aimed to retrieve parameter names, values, and types.
        :param docstring: The docstring to parse.
        :return: A list of tuples containing parameter names, default values, and types.
        """

        parameters = []
        docstring = re.search(r"\((.*?)\)", docstring).group(
            1
        )  # Return class name and init signature
        docstring = docstring.split(",")

        for parameter in docstring:
            if parameter.startswith("self"):
                continue

            name = parameter
            default = None
            parameter_type = "Any"

            if ":" in parameter:
                split_by_semicolon = parameter.split(":", 1)
                name = split_by_semicolon[0].strip()
                type_and_default = split_by_semicolon[1].strip()
                if "=" in type_and_default:
                    split_by_equals = type_and_default.split("=", 1)
                    parameter_type = split_by_equals[0].strip()
                    default = split_by_equals[1].strip()
                    if default.startswith("'") and default.endswith("'"):
                        default = default[1:-1]
                else:
                    parameter_type = type_and_default

            parameters.append((name, default, parameter_type))

        return parameters

    @staticmethod
    def class_parameters_with_defaults(module_name):
        """
        Given a module name, outputs a dictionary of class names and their parameters.
        Keys are class names, and values are lists of parameter information (name, default value, type).
        :param module_name: The module to inspect.
        :return: A dictionary with class names as keys and parameter information as values.
        """

        classes = generalFunctions.find_classes(module_name)
        docstrings = generalFunctions.find_init_docstring_for_classes(classes)

        result = {}

        for class_name, docstring in docstrings.items():
            parameters = generalFunctions.extract_parameters(docstring)
            result[class_name] = parameters

        return result

    @staticmethod
    def select_classes(module_name):
        """
        Given a module name, outputs a list of all class names in the module.
        :param module_name: The module to inspect.
        :return: A list of class names.
        """

        return list(generalFunctions.class_parameters_with_defaults(module_name))

    @staticmethod
    def convert_to_correct_type(value, desired_type):
        """
        Converts the given value to the desired type.
        :param value: The value to convert.
        :param desired_type: The type to convert the value to ('int', 'float', 'str').
        :return: The value converted to the desired type.
        """

        if value is None:
            raise ValueError("Cannot convert to desired type")
        if desired_type == "int":
            return int(value)
        elif desired_type == "float":
            return float(value)
        elif desired_type == "str":
            return str(value)
        else:
            raise ValueError("Unknown type")

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
            elif input_section == "space_charge":
                state.dirty("max_level")

        elif input_section == "all":
            state.update(DashboardDefaults.DEFAULT_VALUES)
            state.dirty("distribution_type")
            state.selected_lattice_list = []
            state.dirty("max_level")
