"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .. import ctrl, state
from .defaults import DashboardDefaults, STATE_INPUTS, determine_section_name
from .utils import GeneralFunctions
from .validation import InputsValidator, errors_tracker

# Set of dropdown input state variables, automatically populated when InputComponents.select() is called
# Used to exclude dropdown inputs from validation since they're constrained to valid options
DROPDOWN_INPUTS = set()


class SharedUtilities:
    @staticmethod
    @state.change(*STATE_INPUTS)
    def on_input_state_change(**_):
        """
        Called when any non-nested state variables are modified.
        """
        non_dropdown_inputs = set(STATE_INPUTS) - DROPDOWN_INPUTS
        state_changes = state.modified_keys & non_dropdown_inputs

        non_dropdown_inputs = set(STATE_INPUTS) - DROPDOWN_INPUTS
        state_changes = state.modified_keys & non_dropdown_inputs

        for state_name in state_changes:
            input = getattr(state, state_name)
            if type(input) is str:
                validation_result = InputsValidator.validate(state_name, input)
                InputsValidator.update_error_message_on_ui(
                    state_name, validation_result
                )

                if not validation_result:
                    GeneralFunctions.set_state_to_numeric(state_name)

                    match state_name:
                        case "kin_energy_on_ui":
                            state.dirty("kin_energy_unit")
                        case _ if (
                            "blocking_factor" in state_name or "n_cell" in state_name
                        ):
                            direction = state_name[-1]
                            InputsValidator.update_n_cell_validation(direction)

                section_name = determine_section_name(state_name)
                errors_tracker.update(section_name, state_name=state_name)

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
