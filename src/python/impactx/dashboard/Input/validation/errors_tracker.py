"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import state
from ..defaults import INPUT_LABELS, DashboardDefaults


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
