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
                if condition == "positive" and value <= 0:
                    errors.append("Must be positive.")
                if condition == "negative" and value >= 0:
                    errors.append("Must be negative.")

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

        n_cell_errors = DashboardValidation.validate_against(n_cell_value, "int")
        blocking_factor_errors = DashboardValidation.validate_against(
            blocking_factor_value, "int", ["non_zero", "positive"]
        )

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
