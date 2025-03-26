"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server
from .. import CardComponents, generalFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


class SpaceChargeFunctions:
    @staticmethod
    def multigrid_settings():
        CardComponents.card_button(
            "mdi-cog",
            color="grey-darken-2",
            click="space_charge_dialog_settings = true",
            v_if="poisson_solver == 'multigrid'",
            documentation="Settings",
        )

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

        n_cell_errors = generalFunctions.validate_against(n_cell_value, "int")
        blocking_factor_errors = generalFunctions.validate_against(
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
