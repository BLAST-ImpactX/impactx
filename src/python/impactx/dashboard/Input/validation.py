"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import keyword
from typing import Union

from .. import state
from .defaults import INPUT_LABELS, DashboardDefaults
from .utils import GeneralFunctions

ALLOWED_INPUT_TYPES = {"int", "float", "str"}
INT_ERROR_MESSAGE = "Must be an integer"
FLOAT_ERROR_MESSAGE = "Must be a float"
NON_ZERO_ERROR = "Must be non-zero"
POSITIVE_ERROR = "Must be positive"
NEGATIVE_ERROR = "Must be negative"
N_CELL_MULTIPLE_ERROR = "Must be a multiple of its blocking factor"
PYTHON_IDENTIFIER_ERROR = "Must be a valid Python identifier"

# Utilized in prob_relative validation
GREATER_THAN_THREE_ERROR = "Must be greater than 3"
GREATER_THAN_ONE_ERROR = "Must be greater than 1"
LESS_THAN_PREVIOUS_ERROR = "Must be less than previous value"


class DashboardValidation:
    """
    Contains all validation logic for the ImpactX dashboard inputs.
    """

    @staticmethod
    def update_error_message_on_ui(state_name: str, error_message: str) -> None:
        """
        Sets the error message for a given state input field.

        :param state_name: The name of the state field to attach the error message to.
        :param error_message: The error message to set.
        """
        validation_name = f"{state_name}_error_message"

        if not hasattr(state, validation_name):
            raise AttributeError(
                f"The error message state does not exist: {validation_name}'"
            )

        if getattr(state, validation_name) != error_message:
            setattr(state, validation_name, error_message)

    @staticmethod
    def validate(
        input_name: str,
        input_value: Union[float, int],
        category: str | None = None,
        parameter_type: str | None = None,
    ) -> list[str]:
        """
        Validates the input value against its default type and any additional conditions.

        :param input_name: The name of the input to validate.
        :param input_value: The value to validate.
        :param category: The category of validation (e.g., 'distribution', 'lattice').
        :param parameter_type: The explicit type to use ('int', 'float', 'str'). If provided, overrides type lookup.
        :return: A list of error messages. An empty list if there are no errors.
        """
        input_type = DashboardValidation._get_input_type(
            input_name, category, parameter_type
        )

        if input_type not in ALLOWED_INPUT_TYPES:
            return [f"Unknown or unsupported type '{input_type}'"]

        if input_type == "str":
            if not DashboardValidation.is_valid_input_name(str(input_value)):
                return [PYTHON_IDENTIFIER_ERROR]
            return []

        numeric_input = GeneralFunctions.convert_to_numeric(input_value)
        type_errors = DashboardValidation._validate_type(numeric_input, input_type)

        if type_errors:
            return type_errors

        additional_validation = DashboardValidation._validate_additional_conditions(
            input_name, numeric_input
        )
        return additional_validation

    @staticmethod
    def _get_input_type(
        input_name: str, category: str | None, parameter_type: str | None
    ) -> str:
        """
        Retrieves the expected input type for validation.

        :param input_name: Name of the parameter.
        :param category: The category of validation (e.g., 'distribution', 'lattice').
        :param parameter_type: The explicit type to use ('int', 'float', 'str'). If provided, overrides type lookup.
        :return: The resolved type as a string ('int', 'float', or 'str').
        """

        if parameter_type:
            return parameter_type

        input_type = None
        if category in {"distribution", "lattice"}:
            input_type = GeneralFunctions.get_default(category, "types")
        else:
            input_type = GeneralFunctions.get_default(input_name, "types")

        return input_type if input_type in ALLOWED_INPUT_TYPES else "str"

    @staticmethod
    def _validate_type(
        numeric_input: Union[float, int, None], value_type: str
    ) -> list[str]:
        """
        Validates a numeric input against the expected type ('int' or 'float').

        :param numeric_input: The value to validate (already converted to numeric or None).
        :param value_type: The expected type ('int' or 'float').
        :return: A list of error messages. Empty if valid.
        """

        if numeric_input is None:
            error_message = (
                INT_ERROR_MESSAGE if value_type == "int" else FLOAT_ERROR_MESSAGE
            )
            return [error_message]

        is_int = isinstance(numeric_input, int)
        is_float = isinstance(numeric_input, (int, float))

        if value_type == "int" and not is_int:
            return [INT_ERROR_MESSAGE]
        elif value_type == "float" and not is_float:
            return [FLOAT_ERROR_MESSAGE]

        return []

    @staticmethod
    def _validate_additional_conditions(
        input_name: str, value: Union[float, int]
    ) -> list[str]:
        """
        Validates additional numeric conditions (e.g., non-zero, positive, negative)
        based on rules defined for the input name.

        :param input_name: Name of the input parameter.
        :param value: Numeric value to validate.
        :return: List of error messages. Empty if all conditions are satisfied.
        """

        lookup_name = "lambda" if "lambda" in input_name else input_name
        additional_conditions = (
            GeneralFunctions.get_default(lookup_name, "validation_condition") or []
        )

        errors = []
        for condition in additional_conditions:
            if condition == "non_zero" and value == 0:
                errors.append(NON_ZERO_ERROR)
            elif condition == "positive" and value <= 0:
                errors.append(POSITIVE_ERROR)
            elif condition == "negative" and value >= 0:
                errors.append(NEGATIVE_ERROR)

        return errors

    @staticmethod
    def is_valid_input_name(user_input: str) -> bool:
        """
        Check if the user input is a valid Python name.
        """
        if user_input is None:
            return True
        return user_input.isidentifier() and not keyword.iskeyword(user_input)

    @staticmethod
    def update_n_cell_validation(direction: str) -> None:
        """
        Validates whether the 'n_cell_<direction>' value is a multiple of 'blocking_factor_<direction>'.

        :param direction: One of 'x', 'y', or 'z' indicating which directional value to check.
        """

        n_cell = GeneralFunctions.convert_to_numeric(
            getattr(state, f"n_cell_{direction}", None)
        )
        blocking_factor = GeneralFunctions.convert_to_numeric(
            getattr(state, f"blocking_factor_{direction}", None)
        )

        if blocking_factor is None or blocking_factor == 0:
            return

        if n_cell % blocking_factor != 0:
            DashboardValidation.update_error_message_on_ui(
                f"n_cell_{direction}", N_CELL_MULTIPLE_ERROR
            )
        else:
            DashboardValidation.update_error_message_on_ui(f"n_cell_{direction}", "")

    @staticmethod
    def validate_prob_relative_fields(index: int, prob_relative_value: float) -> str:
        """
        Validates the prob_relative_fields based on the index and solver type.

        :param index: Index of the modified prob_relative_field.
        :param prob_relative_value: User-provided value to validate.
        :return: Error message if invalid.
        """
        error_message = ""

        try:
            prob_relative_value = float(prob_relative_value)
            poisson_solver = state.poisson_solver

            if index == 0:
                if poisson_solver == "multigrid":
                    if prob_relative_value <= 3:
                        error_message = GREATER_THAN_THREE_ERROR
                elif poisson_solver == "fft":
                    if prob_relative_value <= 1:
                        error_message = GREATER_THAN_ONE_ERROR
            else:
                previous_value = float(state.prob_relative[index - 1])
                if prob_relative_value >= previous_value:
                    error_message = f"{LESS_THAN_PREVIOUS_ERROR} ({previous_value})"
                else:
                    if prob_relative_value <= 1:
                        error_message = GREATER_THAN_ONE_ERROR
        except ValueError:
            error_message = FLOAT_ERROR_MESSAGE

        return error_message


class InputValidator:
    """
    Helps determine whether the simulation can be run based on the current input errors, as well
    as displaying the errors on the UI.
    """

    def __init__(self):
        self.errors = {section: [] for section in DashboardDefaults.INPUT_SECTIONS}

    def _normalize_input_name(self, name: str) -> str:
        """
        Converts a section name to a normalized format suitable for use as a v_model_name.
        EX: "Simulation Parameters" -> "simulation_parameters"

        :param name: The name of a section/input name to normalize.
        """
        return name.lower().replace(" ", "_")

    def _get_error_message(self, input_name: str) -> str:
        """
        Retrieves the error v_model_name for the input_name.

        :param input_name: The name of the input field that is modified.
        """
        normalized = self._normalize_input_name(input_name)
        return getattr(state, f"{normalized}_error_message", "")

    def _get_validation_method(self, section_name: str):
        METHOD_OVERRIDE = {"Simulation Parameters": "_check_input_params_errors"}

        normalized = self._normalize_input_name(section_name)
        method_name = METHOD_OVERRIDE.get(section_name, f"_check_{normalized}_errors")
        return getattr(self, method_name, None), method_name

    def _check_state_errors(self, category: str, input_ui_name: str) -> None:
        """
        Appends any error for a UI input to its category.

        :param category: The dashboard section name.
        :param input_ui_name: The UI label (e.g., "CSR Bins", "Max Level").
        """
        error = self._get_error_message(input_ui_name)
        print(error)
        if error:
            self.errors[category].append(f"{input_ui_name}: {error}")

    def _check_input_params_errors(self) -> list[str]:
        errors: list[str] = []
        for v_model_name, display_name in INPUT_LABELS.items():
            if self._get_error_message(v_model_name):
                errors.append(
                    f"{display_name}: {self._get_error_message(v_model_name)}"
                )
        return errors

    def _check_space_charge_errors(self) -> list[str]:
        # If space charge is disabled, return no errors
        if getattr(state, "space_charge", None) == "false":
            return []

        errors: list[str] = []
        directions = ["x", "y", "z"]
        inputs = ["n_cell", "blocking_factor"]

        for axis in directions:
            for param in inputs:
                name = f"{param}_{axis}"
                error = self._get_error_message(name)
                if error:
                    errors.append(f"{name}: {error}")

        for i, field in enumerate(state.prob_relative_fields):
            if field.get("error_message"):
                errors.append(f"prob_relative[{i}]: {field['error_message']}")

        return errors

    def _check_variables_errors(self) -> list[str]:
        for variables in state.variables:
            if variables.get("error_message"):
                return ["Variable definition error"]
        return []

    def _check_csr_errors(self) -> list[str]:
        # Clear previous CSR errors
        self.errors["CSR"] = []

        # If CSR is disabled, return no errors
        if not getattr(state, "csr", False):
            return []

        self._check_state_errors("CSR", "csr_bins")
        return self.errors["CSR"]

    def _check_distribution_parameters_errors(self) -> list[str]:
        errors = []
        for name, param in state.selected_distribution_parameters.items():
            if param["error_message"]:
                errors.append(f"{name}: {param['error_message']}")
        return errors

    def _check_lattice_configuration_errors(self) -> list[str]:
        errors: list[str] = []
        for index, lattice in enumerate(state.selected_lattice_list):
            for param in lattice["parameters"]:
                if param["parameter_error_message"]:
                    errors.append(
                        f"Element {index + 1} ({lattice['name']}) - {param['parameter_name']}: {param['parameter_error_message']}"
                    )
        if not state.selected_lattice_list:
            errors.append("Lattice is empty")

        periods_error = getattr(state, "periods_error_message", "")
        if periods_error:
            errors.append(f"Periods: {periods_error}")

        return errors

    def update(self, section_name: str) -> None:
        """
        Updates the state that contains all input errors displaying on the dashboard.

        Helpful for determining if a simulation can be ran on the dashboard and if/what are the input errors.
        """
        section_errors = []
        checker_method, expected_method_name = self._get_validation_method(section_name)

        if callable(checker_method):
            section_errors = checker_method()
        else:
            print(
                f"[WARNING] No validation method found for section '{section_name}'. "
                f"Expected method: '{expected_method_name}'"
            )

        self.errors[section_name] = section_errors

        categorized_errors = [
            {"category": category, "errors": errors}
            for category, errors in self.errors.items()
            if errors
        ]

        state.number_of_input_errors = sum(
            len(item["errors"]) for item in categorized_errors
        )

        state.disable_simulation = bool(categorized_errors)
        state.input_errors = categorized_errors


sim_validator = InputValidator()
