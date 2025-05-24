"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .. import state

from . import generalFunctions

ALLOWED_INPUT_TYPES = {"int", "float", "str"}
INT_ERROR_MESSAGE = "Must be an integer"
FLOAT_ERROR_MESSAGE = "Must be a float"
NON_ZERO_ERROR = "Must be non-zero."
POSITIVE_ERROR = "Must be positive."
NEGATIVE_ERROR = "Must be negative."

class DashboardValidation:
    """
    Contains all validation logic for the ImpactX dashboard inputs.
    """

    @staticmethod
    def validate_input(input_name: str, input_value, category=None, parameter_type=None):
        """
        Validates the input value against the desired type and additional conditions.

        :param input_name: The name of the parameter being validated.
        :param input_value: The value to validate.
        :param category: The category of validation (e.g., 'distribution', 'lattice').
        :param parameter_type: The explicit type to use ('int', 'float', 'str'). If provided, overrides type lookup.
        :return: A list of error messages. An empty list if there are no errors.
        """
        input_type = DashboardValidation._get_input_type(input_name, category, parameter_type)

        if input_type not in ALLOWED_INPUT_TYPES:
            return [f"Unknown or unsupported type '{input_type}'"]

        if input_type == "str":
            no_errors = []
            return no_errors
        
        numeric_input = generalFunctions.convert_to_numeric(input_value)
        value, type_errors = DashboardValidation._validate_type(numeric_input, input_type)

        if type_errors:
            return type_errors

        additional_validation = DashboardValidation._validate_additional_conditions(input_name, value)
        return additional_validation

    @staticmethod
    def _get_input_type(input_name: str, category, parameter_type):
        """
        Helper method to determine the input type.
        """
        if parameter_type is not None:
            return parameter_type
    
        if category in ["distribution", "lattice"]:
            input_type = "float"
        else:
            input_type = generalFunctions.get_default(input_name, "types")

        return input_type

    @staticmethod
    def _validate_type(numeric_input, value_type):
        """
        Helper method to validate the numeric_input.
        """
        if numeric_input is None:
            error_message = (
                INT_ERROR_MESSAGE if value_type == "int" else FLOAT_ERROR_MESSAGE
            )
            return None, [error_message]

        is_int = isinstance(numeric_input, int)
        is_float = isinstance(numeric_input, (int, float))

        if value_type == "int" and not is_int:
            return None, [INT_ERROR_MESSAGE]
        elif value_type == "float" and not is_float:
            return None, [FLOAT_ERROR_MESSAGE]

        return numeric_input, []

    @staticmethod
    def _validate_additional_conditions(input_name: str, value):
        """
        Helper method to validate additional conditions (ie. non-zero, positive, negative).
        """

        if value is None:
            return []
            
        lookup_name = "lambda" if "lambda" in input_name else input_name
        additional_conditions = generalFunctions.get_default(lookup_name, "validation_condition") or []
        
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
    def validate_prob_relative_fields(index, prob_relative_value):
        """
        This function checks specific validation requirements
        for prob_relative_fields.
        :param index: The index of the prob_relative_field modified.
        :param prob_relative_value: The numerical value entered by the user.
        :return: An error message. An empty string if there is no error.
        """
        error_message = ""

        try:
            prob_relative_value = float(prob_relative_value)
            poisson_solver = state.poisson_solver

            if index == 0:
                if poisson_solver == "multigrid":
                    if prob_relative_value <= 3:
                        error_message = "Must be greater than 3."
                elif poisson_solver == "fft":
                    if prob_relative_value <= 1:
                        error_message = "Must be greater than 1."
            else:
                previous_value = float(state.prob_relative[index - 1])
                if prob_relative_value >= previous_value:
                    error_message = (
                        f"Must be less than previous value ({previous_value})."
                    )
                else:
                    if prob_relative_value <= 1:
                        error_message = "Must be greater than 1."
        except ValueError:
            error_message = "Must be a float."

        return error_message

    @staticmethod
    def validate_n_cell_and_blocking_factor(direction):
        """
        Validation function for n_cell and blocking_factor parameters.
        """
        n_cell_value = getattr(state, f"n_cell_{direction}", None)
        blocking_factor_value = getattr(state, f"blocking_factor_{direction}", None)

        n_cell_errors = DashboardValidation.validate_input("n_cell", n_cell_value)
        blocking_factor_errors = DashboardValidation.validate_input("blocking_factor", blocking_factor_value)

        setattr(state, f"n_cell_{direction}_error_message", "; ".join(n_cell_errors))
        setattr(
            state,
            f"blocking_factor_{direction}_error_message",
            "; ".join(blocking_factor_errors),
        )

        if not n_cell_errors and not blocking_factor_errors:
            n_cell_value = int(n_cell_value)
            blocking_factor_value = int(blocking_factor_value)
            if n_cell_value % blocking_factor_value != 0:
                setattr(
                    state,
                    f"n_cell_{direction}_error_message",
                    "Must be a multiple of blocking factor.",
                )

    @staticmethod
    def update_simulation_validation_status():
        """
        Checks if any input fields are not provided with the correct input type.
        Updates the state to enable or disable the run simulation button.
        """

        error_details = []

        # Check for errors in distribution parameters
        for param_name, param in state.selected_distribution_parameters.items():
            if param["error_message"]:
                error_details.append(f"{param_name}: {param['error_message']}")

        # Check for errors in lattice parameters
        for lattice in state.selected_lattice_list:
            for param in lattice["parameters"]:
                if param["parameter_error_message"]:
                    error_details.append(
                        f"Lattice {lattice['name']} - {param['parameter_name']}: {param['parameter_error_message']}"
                    )

        # Check for errors in input card
        if state.npart_error_message:
            error_details.append(f"Number of Particles: {state.npart_error_message}")
        if state.kin_energy_error_message:
            error_details.append(f"Kinetic Energy: {state.kin_energy_error_message}")
        if state.bunch_charge_C_error_message:
            error_details.append(f"Bunch Charge: {state.bunch_charge_C_error_message}")
        if state.charge_qe_error_message:
            error_details.append(
                f"Ref. Particle Charge: {state.charge_qe_error_message}"
            )
        if state.mass_MeV_error_message:
            error_details.append(f"Ref. Particle Mass: {state.mass_MeV_error_message}")

        if state.selected_lattice_list == []:
            error_details.append("LatticeListIsEmpty")
        if state.periods_error_message:
            error_details.append(f"Periods: {state.periods_error_message}")

        # Check for errors in CSR parameters
        if state.csr_bins_error_message:
            error_details.append(f"CSR Bins: {state.csr_bins_error_message}")

        # Check for errors in Space Charge parameters
        if state.space_charge:
            # n_cell parameters
            for direction in ["x", "y", "z"]:
                n_cell_error = getattr(state, f"error_message_n_cell_{direction}")
                if n_cell_error:
                    error_details.append(f"n_cell_{direction}: {n_cell_error}")

            # Blocking factor parameters
            for direction in ["x", "y", "z"]:
                blocking_factor_error = getattr(
                    state, f"error_message_blocking_factor_{direction}"
                )
                if blocking_factor_error:
                    error_details.append(
                        f"blocking_factor_{direction}: {blocking_factor_error}"
                    )

            # Prob Relative Fields
            for index, field in enumerate(state.prob_relative_fields):
                if field["error_message"]:
                    error_details.append(
                        f"prob_relative[{index}]: {field['error_message']}"
                    )

        def has_error_in_variables() -> bool:
            """
            Determines if state.variables contains an error message.
            Return true if yes, false if no. Needed to not allow sim. to run
            if there is an error.
            """
            results = any(
                variable.get("error_message", "") for variable in state.variables
            )
            return results

        if has_error_in_variables():
            error_details.append("error")

        state.disableRunSimulationButton = bool(error_details)
