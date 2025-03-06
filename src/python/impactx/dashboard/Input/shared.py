"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .. import setup_server
from ..Input.inputParameters.inputMain import InputParameters
from . import DashboardDefaults, generalFunctions

server, state, ctrl = setup_server()


input_parameters_defaults = list(DashboardDefaults.INPUT_PARAMETERS.keys())
space_charge_defaults = list(DashboardDefaults.CSR.keys())
INPUT_DEFAULTS = input_parameters_defaults + space_charge_defaults


class SharedUtilities:
    @staticmethod
    @state.change(*INPUT_DEFAULTS)
    def on_input_state_change(**_):
        state_changes = state.modified_keys & set(INPUT_DEFAULTS)
        for state_name in state_changes:
            if type(state[state_name]) is str:
                value = getattr(state, state_name)
                desired_type = DashboardDefaults.TYPES.get(state_name, None)
                validation_name = f"{state_name}_error_message"
                conditions = DashboardDefaults.VALIDATION_CONDITION.get(
                    state_name, None
                )

                validation_result = generalFunctions.validate_against(
                    value, desired_type, conditions
                )
                setattr(state, validation_name, validation_result)
                generalFunctions.update_simulation_validation_status()

                if validation_result == []:
                    converted_value = generalFunctions.convert_to_correct_type(
                        value, desired_type
                    )

                    if getattr(state, state_name) != converted_value:
                        setattr(state, state_name, converted_value)
                        if state_name == "kin_energy_on_ui":
                            InputParameters.on_kin_energy_unit_change()

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
