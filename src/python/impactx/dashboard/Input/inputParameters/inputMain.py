"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardComponents, DashboardDefaults, InputComponents, generalFunctions
from . import InputFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Callbacks
# -----------------------------------------------------------------------------


input_parameters_defaults = list(DashboardDefaults.INPUT_PARAMETERS.keys())
space_charge_defaults = list(DashboardDefaults.CSR.keys())
INPUT_DEFAULTS = input_parameters_defaults + space_charge_defaults


@state.change(*INPUT_DEFAULTS)
def on_input_state_change(**_):
    state_changes = state.modified_keys & set(INPUT_DEFAULTS)
    for state_name in state_changes:
        if type(state[state_name]) is str:
            print(
                f"{state_name} = {state[state_name]} and type: {type(state[state_name])}"
            )
            value = getattr(state, state_name)
            desired_type = DashboardDefaults.TYPES.get(state_name, None)
            validation_name = f"{state_name}_error_message"
            conditions = DashboardDefaults.VALIDATION_CONDITION.get(state_name, None)

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
                        on_kin_energy_unit_change()


@state.change("kin_energy_unit")
def on_kin_energy_unit_change(**kwargs) -> None:
    if state.kin_energy_on_ui != 0:
        InputFunctions.update_kin_energy_sim_value()


# -----------------------------------------------------------------------------
# Content
# -----------------------------------------------------------------------------


class InputParameters:
    """
    User-Input section for beam properties.
    """

    def card(self):
        """
        Creates UI content for beam properties.
        """

        with vuetify.VCard(style="width: 340px; height: 350px"):
            CardComponents.input_header("Input Parameters")
            with vuetify.VCardText():
                with vuetify.VRow(classes="py-2"):
                    with vuetify.VCol(cols=6, classes="py-0"):
                        vuetify.VCheckbox(
                            label="Space Charge",
                            v_model=("space_charge", False),
                            dense=True,
                        )
                    with vuetify.VCol(cols=6, classes="py-0"):
                        vuetify.VCheckbox(
                            label="CSR",
                            v_model=("csr", False),
                            dense=True,
                        )
                with vuetify.VRow(classes="my-2"):
                    with vuetify.VCol(cols=6, classes="py-0"):
                        InputComponents.text_field(
                            label="Ref. Particle Charge",
                            v_model_name="charge_qe",
                        )
                    with vuetify.VCol(cols=6, classes="py-0"):
                        InputComponents.text_field(
                            label="Ref. Particle Mass",
                            v_model_name="mass_MeV",
                        )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=12, classes="py-0"):
                        InputComponents.text_field(
                            label="Number of Particles",
                            v_model_name="npart",
                        )
                with vuetify.VRow(classes="my-2"):
                    with vuetify.VCol(cols=8, classes="py-0"):
                        InputComponents.text_field(
                            label="Kinetic Energy",
                            v_model_name="kin_energy_on_ui",
                            classes="mr-2",
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        InputComponents.select(
                            label="Unit",
                            v_model_name="kin_energy_unit",
                        )
                with vuetify.VRow(classes="my-2"):
                    with vuetify.VCol(cols=12, classes="py-0"):
                        InputComponents.text_field(
                            label="Bunch Charge",
                            v_model_name="bunch_charge_C",
                        )
