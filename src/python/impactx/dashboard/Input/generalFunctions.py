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
    def convert_to_numeric(input: str) -> int | float:
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

            if "Optional" in parameter_type:
                parameter_type = parameter_type[len("Optional[") : -1]
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
