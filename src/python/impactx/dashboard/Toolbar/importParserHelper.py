"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import ast
import re

from ..Input.defaults import DashboardDefaults


class DashboardParserHelper:
    """
    Helper functions for building dashboard parser.
    """

    @staticmethod
    def import_file_content(file: str, state) -> dict:
        """
        Retrieves and prints the content of the uploaded simulation file.

        :param content: The content of the ImpactX simulation file.
        """
        if file:
            content = file["content"].decode("utf-8")
            return content
        else:
            state.file_content = ""
            return ""

    @staticmethod
    def parse_single_inputs(content: str) -> dict:
        """
        Parses individual input parameters from the simulation file content.

        :param content: The content of the ImpactX simulation file.
        """
        reference_dictionary = DashboardDefaults.DEFAULT_VALUES.copy()

        parsing_patterns = [
            r"\b{}\s*=\s*([^#\n]+)",  # (param = value)
            r"set_{}\(([^)]+)\)",  # (set_param(value))
        ]

        for parameter_name in reference_dictionary.keys():
            if parameter_name.endswith("_list"):
                continue

            for pattern in parsing_patterns:
                pattern_match = re.search(pattern.format(parameter_name), content)
                if pattern_match:
                    value = ast.literal_eval(pattern_match.group(1))
                    reference_dictionary[parameter_name] = value
                    break

        # Handling for kin_energy
        kin_energy_pattern_match = re.search(
            r"\bkin_energy_MeV\s*=\s*([^#\n]+)", content
        )
        if kin_energy_pattern_match:
            kin_energy_value = kin_energy_pattern_match.group(1)
            reference_dictionary["kin_energy_on_ui"] = kin_energy_value

        return reference_dictionary

    @staticmethod
    def parse_list_inputs(content: str) -> dict:
        """
        Parses list-based input parameters from the simulation file content.

        :param content: The content of the ImpactX simulation file.
        """
        dictionary = {}
        list_inputs = ["n_cell", "prob_relative"]
        list_parsing = "{} = (\\[.*?\\])"

        for input_name in list_inputs:
            match = re.search(list_parsing.format(input_name), content)
            if match:
                values = ast.literal_eval(match.group(1).strip())

                if input_name == "n_cell":
                    for i, dim in enumerate(["x", "y", "z"]):
                        dictionary[f"n_cell_{dim}"] = values[i]

                if input_name == "prob_relative":
                    dictionary["prob_relative"] = values

        return dictionary

    @staticmethod
    def parse_distribution(content: str) -> dict:
        """
        Parses distribution section from the simulation file content.

        :param content: The content of the ImpactX simulation file.
        """

        dictionary = {"distribution": {"name": "", "type": "", "parameters": {}}}

        distribution_name = re.search(r"distribution\.(\w+)\(", content)
        distribution_type_twiss = re.search(r"twiss\((.*?)\)", content, re.DOTALL)
        distribution_type_quadratic = re.search(
            r"distribution\.\w+\((.*?)\)", content, re.DOTALL
        )
        parameters = {}

        def extract_parameters(distribution_type, parsing_pattern):
            parameter_pairs = re.findall(parsing_pattern, distribution_type.group(1))
            parsed_parameters = {}

            for param_name, param_value in parameter_pairs:
                parsed_parameters[param_name] = param_value
            return parsed_parameters

        if distribution_name:
            dictionary["distribution"]["name"] = distribution_name.group(1)

        if distribution_type_twiss:
            dictionary["distribution"]["type"] = "Twiss"
            parameters = extract_parameters(
                distribution_type_twiss, r"(\w+)=(\d+\.?\d*)"
            )
        elif distribution_type_quadratic:
            dictionary["distribution"]["type"] = "Quadratic"
            parameters = extract_parameters(
                distribution_type_quadratic, r"(\w+)=([^,\)]+)"
            )

        dictionary["distribution"]["parameters"] = parameters

        return dictionary

    @staticmethod
    def parse_lattice_elements(content: str) -> dict:
        """
        Parses lattice elements from the simulation file content.

        :param content: The content of the ImpactX simulation file.
        """

        dictionary = {"lattice_elements": []}

        lattice_elements = re.findall(r"elements\.(\w+)\((.*?)\)", content)

        for element_name, element_parameter in lattice_elements:
            element = {"element": element_name, "parameters": {}}

            parameter_pairs = re.findall(r"(\w+)=([^,\)]+)", element_parameter)
            for parameter_name, parameter_value in parameter_pairs:
                parameter_value_cleaned = parameter_value.strip("'\"")
                element["parameters"][parameter_name] = parameter_value_cleaned

            dictionary["lattice_elements"].append(element)

        return dictionary
