"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import inspect
from typing import Any, List, Tuple

from impactx.distribution_input_helpers import twiss

from ... import state
from .. import GeneralFunctions


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
            param_name: float(param["value"]) if param_is_valid else 0.0
            for param_name, param in state.selected_distribution_parameters.items()
            if (param_is_valid := param["error_message"] == [])
        }

        return parameter_input

    @staticmethod
    def get_distribution_units(name: str) -> str:
        """
        Returns the correct units depending on if
        selected_distribution == Twiss.
        """
        if "beta" in name or "emitt" in name:
            return GeneralFunctions.get_default(name, "units")
        return ""

    @staticmethod
    def get_twiss_data() -> List[Tuple[str, Any, type]]:
        """
        Retrieves parameters names and default values for the Twiss parameters.

        Utilizes the twiss helper function from `distribution_input_helpers`.
        """
        param_data = []

        sig = inspect.signature(twiss)
        for parameter in sig.parameters.values():
            name = parameter.name
            default_value = (
                parameter.default if parameter.default != inspect._empty else None
            )
            default_type = GeneralFunctions.get_default(name, "types")
            param_data.append((name, default_value, default_type))
        return param_data
