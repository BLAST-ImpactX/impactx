"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import ast
import re

from ....Input.defaults import DashboardDefaults


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
        Parses individual simulation parameters from the simulation file content.

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

                    if parameter_name == "space_charge" and isinstance(value, bool):
                        value = str(value).lower()

                    reference_dictionary[parameter_name] = value
                    break

        # Handling for kin_energy
        kin_energy_pattern_match = re.search(
            r"\bkin_energy_MeV\s*=\s*([^#\n]+)", content
        )
        if kin_energy_pattern_match:
            kin_energy_value = kin_energy_pattern_match.group(1).strip()

            try:
                parsed_value = ast.literal_eval(kin_energy_value)
            except (ValueError, SyntaxError):
                parsed_value = reference_dictionary.get("kin_energy_MeV")

            reference_dictionary["kin_energy_on_ui"] = parsed_value

        return reference_dictionary

    @staticmethod
    def parse_list_inputs(content: str) -> dict:
        """
        Parses list-based simulation parameters from the simulation file content.

        Some of the dashboard inputs are lists, such as n_cell = [16, 16, 16].
        We have to look inside of the brackets to retrieve the individual values.

        :param content: The content of the ImpactX simulation file.
        """
        dictionary = {}
        list_inputs = [
            "n_cell",
            "prob_relative",
            "blocking_factor_x",
            "blocking_factor_y",
            "blocking_factor_z",
        ]

        # The below regex will capture everything inside of the brackets
        # it assumes that the list will be placed in a single line (not multiple).
        # Example: n_cell = [16, 16, 16] -->  [16, 16, 16]
        list_parsing = r"{}\s*=\s*(\[[^\]]*\])"

        for input_name in list_inputs:
            match = re.search(list_parsing.format(input_name), content)
            if match:
                values = ast.literal_eval(match.group(1).strip())

                if input_name == "n_cell":
                    for i, dim in enumerate(["x", "y", "z"]):
                        dictionary[f"n_cell_{dim}"] = values[i]

                elif input_name == "prob_relative":
                    dictionary["prob_relative"] = values

                elif input_name.startswith("blocking_factor_"):
                    # Convert [16] -> 16
                    dictionary[input_name] = values[0]

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
    def parse_variables(content: str, used_vars: set) -> dict:
        """
        Parses variable definitions from the simulation file content.

        The function looks for lines that define variables in similar format of:
            variable_name = value

        variables_regex breakdown:
            ^\s* - Remove leading whitespace before variable name
            (\w+) => Retrieve the variable name
            \s*=\s*  => Allow the line to have spaces around '=' (ie. 'a= 10' or 'a =10')
            ([^#\n]+) => Captures everything after '=' up to a comment '#' or end of line

        EX when calling this function:
            content = '''
                a = 10
                b = 2.5  # some comment
                c = "some string"
            '''
            used_vars = {"a", "b"}
            result = parse_variables(content, used_vars)
            # result: { "a": 10, "b": 2.5}

        :param content: The content of the ImpactX simulation file.
        :param used_vars: Set of variable names to filter by.
        :return: Dictionary of parsed variables.
        """

        variables = {}
        variables_regex = r"^\s*(\w+)\s*=\s*([^#\n]+)"
        variables_found = re.findall(variables_regex, content, re.MULTILINE)

        for var_name, var_value in variables_found:
            if var_name not in used_vars:
                continue
            try:
                value = ast.literal_eval(var_value.strip())
                # only allowing numbers at the moment
                if isinstance(value, (int, float)):
                    variables[var_name] = value
            except Exception:
                continue
        return variables
