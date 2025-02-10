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
