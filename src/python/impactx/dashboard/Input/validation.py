"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .. import setup_server

server, state, ctrl = setup_server()


class DashboardValidation:
    """
    Contains all validation logic for the ImpactX dashboard inputs.
    """

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
