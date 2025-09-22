"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import state
from ..defaults import INPUT_LABELS, DashboardDefaults, determine_section_name, csr_defaults
from ..utils import GeneralFunctions

class ErrorsTracker:
    """
    Helps track the input errors across the dashboard.
    """

    def __init__(self):
        self.errors = {section: [] for section in DashboardDefaults.INPUT_SECTIONS}

    def _get_error_message(self, input_name: str) -> str:
        """
        Retrieves the error v_model_name for the input_name.

        :param input_name: The name of the input field that is modified.
        """
        candidate_names = [input_name]

        normalized = GeneralFunctions.normalize_for_v_model(input_name)
        if normalized not in candidate_names:
            candidate_names.append(normalized)

        for candidate in candidate_names:
            state_attr = f"{candidate}_error_message"
            if hasattr(state, state_attr):
                return getattr(state, state_attr)

        return ""

    def clear_category(self, category: str) -> None:
        self.errors[category] = []
        self._update_ui()

    def _clear_error(self, category: str, input_name: str) -> None:
        """
        Clears any error messages for a specific input in the errors tracker.

        :param category: The dashboard section name.
        :param input_name: The UI label (e.g., "CSR Bins", "Max Level").
        """
        self.errors[category] = [
            msg for msg in self.errors[category] if not msg.startswith(f"{input_name}:")
        ]

    def _get_validation_method(self, section_name: str):
        """
        Returns the validation method for a section and the expected method name.
        """
        method_override = {
            "Simulation Parameters": "_check_input_params_errors",
        }

        normalized = GeneralFunctions.normalize_for_v_model(section_name)
        method_name = method_override.get(section_name, f"_check_{normalized}_errors")

        return getattr(self, method_name, None), method_name

    def _check_state_errors(self, input_name: str) -> None:
        """
        Appends any error for a UI input to its category.

        :param category: The dashboard section name.
        :param input_name: The UI label (e.g., "CSR Bins", "Max Level").
        """
        category = determine_section_name(input_name)
        self._clear_error(category, input_name)
        error = self._get_error_message(input_name)
        if error:
            self.errors[category].append(f"{input_name}: {error}")

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

        # Check MLMG settings if using multigrid poisson solver
        if getattr(state, "poisson_solver", None) == "multigrid":
            mlmg_fields = [
                "mlmg_relative_tolerance",
                "mlmg_absolute_tolerance",
                "mlmg_max_iters",
                "mlmg_verbosity",
            ]
            for field in mlmg_fields:
                error = self._get_error_message(field)
                if error:
                    errors.append(f"{field}: {error}")

        return errors

    def _check_variables_errors(self) -> list[str]:
        for variables in state.variables:
            if variables.get("error_message"):
                return ["Variable definition error"]
        return []

    def _check_csr_errors(self) -> list[str]:
        """
        Returns a list of current CSR-related input errors.
        Only called when the csr section is enabled/disabled or a csr input is modified
        """
        if getattr(state, "csr", None) in ("false", False):
            return []

        errors = []
        for field_name in csr_defaults:
            error_message = self._get_error_message(field_name)
            if isinstance(error_message, list):
                error_message = "; ".join(error_message)
            if error_message:
                errors.append(f"{field_name}: {error_message}")
        return errors

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

    def update(self, section_name: str, state_name: str | bool = False) -> None:
        checker_method, expected_method_name = self._get_validation_method(section_name)

        if isinstance(state_name, str) and state_name:
            self._check_state_errors(state_name)
        else:
            if callable(checker_method):
                self.errors[section_name] = checker_method()
            else:
                print(
                    f"[WARNING] No validation method found for section '{section_name}'. "
                    f"Expected method: '{expected_method_name}'"
                )
                self.errors.setdefault(section_name, [])

        self._update_ui()

    def _update_ui(self) -> None:
        """
        Rebuilds categorized error data and updates the dashboard state counters/flags.
        """
        categorized_errors = []
        total_errors = 0

        for category, errors in self.errors.items():
            if errors:
                categorized_errors.append({"category": category, "errors": errors})
                total_errors += len(errors)

        state.number_of_input_errors = total_errors
        state.disable_simulation = bool(categorized_errors)
        state.input_errors_list = categorized_errors


errors_tracker = ErrorsTracker()
