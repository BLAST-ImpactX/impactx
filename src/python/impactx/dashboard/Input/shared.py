"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .. import ctrl, state
from . import DashboardDefaults, DashboardValidation
from .utils import GeneralFunctions

simulation_parameters_defaults = list(DashboardDefaults.SIMULATION_PARAMETERS.keys())
csr_defaults = list(DashboardDefaults.CSR.keys())
sc_multigrid_defaults = [
    "mlmg_relative_tolerance",
    "mlmg_absolute_tolerance",
    "mlmg_max_iters",
    "mlmg_verbosity",
]

lattice_state_defaults = ["periods"]
INPUT_DEFAULTS = (
    csr_defaults
    + simulation_parameters_defaults
    + sc_multigrid_defaults
    + lattice_state_defaults
)


def update_error_message_on_ui(state_name: str, error_message: str) -> None:
    """
    Called when we want to set an error message to an input
    """
    validation_name = f"{state_name}_error_message"
    setattr(state, validation_name, error_message)

def set_input_to_numeric(state_name: str) -> None:
    """
    Retrieves the value of state_name, converts to a numeric and re-sets the state_name.
    """
    current_input = getattr(state, state_name)
    numeric_input = generalFunctions.convert_to_numeric(current_input)
    setattr(state, state_name, numeric_input)

class SharedUtilities:
    @staticmethod
    @state.change(*INPUT_DEFAULTS)
    def on_input_state_change(**_):
        """
        Called when any non-nested state variables are modified.
        """
        state_changes = state.modified_keys & set(INPUT_DEFAULTS)
        for state_name in state_changes:
            if type(state[state_name]) is str:
                input = getattr(state, state_name)

                validation_result = DashboardValidation.validate_input(
                    state_name, input
                )
                update_error_message_on_ui(state_name, validation_result)

                if not validation_result:
                    set_input_to_numeric(state_name)
                    if state_name == "kin_energy_on_ui":
                        SimulationParameters.on_kin_energy_unit_change()

                DashboardValidation.update_simulation_validation_status()

    @ctrl.add("collapse_all_sections")
    def on_collapse_all_sections_click():
        state.expand_all_sections = not state.expand_all_sections
        for collapsable_section in DashboardDefaults.COLLAPSABLE_SECTIONS:
            setattr(state, collapsable_section, state.expand_all_sections)

    @state.change(*DashboardDefaults.COLLAPSABLE_SECTIONS)
    def on_collapsable_section_change(**kwargs):
        max_height = "1000px"
        min_height = "3.75rem"

        state_changes = state.modified_keys & set(
            DashboardDefaults.COLLAPSABLE_SECTIONS
        )
        for state_name in state_changes:
            new_height = min_height if getattr(state, state_name) else max_height

            setattr(
                state,
                f"{state_name}_height",
                {
                    "max-height": new_height,
                    "overflow": "hidden",
                    "transition": "max-height 0.5s",
                },
            )
