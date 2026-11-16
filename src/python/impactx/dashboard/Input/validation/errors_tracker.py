"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import state
from ..defaults import beam_intensity_name


class ErrorsTracker:
    def update_simulation_validation_status(self):
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
        kin_energy_error = getattr(state, "kin_energy_on_ui_error_message", None)
        if kin_energy_error:
            error_details.append(f"Kinetic Energy: {kin_energy_error}")
        # only the beam intensity input shown for the tracking mode can block a run
        if beam_intensity_name() == "beam_current_A":
            if state.beam_current_A_error_message:
                error_details.append(
                    f"Beam Current: {state.beam_current_A_error_message}"
                )
        elif state.bunch_charge_C_error_message:
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
        if state.csr and state.csr_bins_error_message:
            error_details.append(f"CSR Bins: {state.csr_bins_error_message}")

        # Check for errors in Space Charge parameters
        if state.space_charge != "false":
            # grid parameters
            for field in ["n_cell", "blocking_factor"]:
                for direction in ["x", "y", "z"]:
                    name = f"{field}_{direction}"
                    error = getattr(state, f"{name}_error_message", None)
                    if error:
                        error_details.append(f"{name}: {error}")

            # Multigrid parameters
            if state.poisson_solver == "multigrid":
                for name in [
                    "mlmg_relative_tolerance",
                    "mlmg_absolute_tolerance",
                    "mlmg_max_iters",
                    "mlmg_verbosity",
                ]:
                    error = getattr(state, f"{name}_error_message", None)
                    if error:
                        error_details.append(f"{name}: {error}")

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


errors_tracker = ErrorsTracker()
