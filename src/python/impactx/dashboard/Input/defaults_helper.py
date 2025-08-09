"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import inspect
import re
from typing import Callable, Dict, List, Type

from impactx.distribution_input_helpers import twiss


class InputDefaultsHelper:
    """
    Methods in this class are used to dynamically parse
    core ImpactX data (default values, docstrings, etc.)
    """

    @staticmethod
    def get_docstrings(
        class_names: List[Type], default_list: Dict[str, any]
    ) -> Dict[str, str]:
        """
        Retrieves docstrings for each method and property
        in the provided clases.

        :param classes: The class names to parse docstrings with.
        :param defaults_list: The dictionary of defaults value.
        """

        docstrings = {}

        for each_class in class_names:
            for name, attribute in inspect.getmembers(each_class):
                if name not in default_list:
                    continue

                is_method = inspect.isfunction(attribute)
                is_property = inspect.isdatadescriptor(attribute)

                if is_method or is_property:
                    docstring = inspect.getdoc(attribute) or ""
                    docstrings[name] = docstring

        distribution_tooltips = InputDefaultsHelper.get_tooltips_from_param(twiss)
        docstrings.update(distribution_tooltips)

        return docstrings

    @staticmethod
    def get_tooltips_from_param(function: Callable) -> Dict[str, str]:
        """
        Extract all ':param name: description' entries from a function's docstring.

        Example:
            :param beta_x: Beta function value in the x dimension.
            :param emitt_x: Emittance function value in the x dimension.

        This will produce:
            {
                "beta_x": "Beta function value in the x dimension.",
                "emitt_x": "Emittance function value in the x dimension."
            }

        :param function: The function whose docstring you want to parse.
        :return: A dict mapping each parameter name to its description.
        """
        tooltip_results = {}
        docstring = inspect.getdoc(function) or ""
        pattern = re.compile(r"^\s*:param\s+(\w+)\s*:\s*(.+)$", re.MULTILINE)
        pattern_matches = list(pattern.finditer(docstring))

        if not pattern_matches:
            raise ValueError(
                f"Found no docstrings to parse in function {function.__name__}"
            )

        for match in pattern_matches:
            param_name = match.group(1)
            param_description = match.group(2)
            tooltip_results[param_name] = param_description

        return tooltip_results

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

            match parameter_type:
                case optional_type if "Optional" in optional_type:
                    parameter_type = parameter_type[len("Optional[") : -1]
                case "typing.SupportsFloat":
                    parameter_type = "float"
                case "typing.SupportsInt":
                    parameter_type = "int"
                case "str | None":
                    parameter_type = "str"

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

        classes = InputDefaultsHelper.find_classes(module_name)
        docstrings = InputDefaultsHelper.find_init_docstring_for_classes(classes)

        result = {}

        for class_name, docstring in docstrings.items():
            parameters = InputDefaultsHelper.extract_parameters(docstring)
            result[class_name] = parameters

        return result

    @staticmethod
    def select_classes(module_name):
        """
        Given a module name, outputs a list of all class names in the module.
        :param module_name: The module to inspect.
        :return: A list of class names.
        """

        return list(InputDefaultsHelper.class_parameters_with_defaults(module_name))
