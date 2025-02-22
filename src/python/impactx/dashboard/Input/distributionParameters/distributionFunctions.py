"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


class DistributionFunctions:
    """
    Helper functions for the distribution parameters.
    """

    @staticmethod
    def convert_distribution_parameters_to_valid_type():
        """
        Helper function to convert user-inputted distribution parameters
        from string type to float type.

        :return: A dictionary with parameter names as keys and their validated values.
        """

        parameter_input = {
            param["parameter_name"]: float(param["parameter_default_value"])
            if param_is_valid
            else 0.0
            for param in state.selected_distribution_parameters
            if (param_is_valid := param["parameter_error_message"] == [])
        }

        return parameter_input
