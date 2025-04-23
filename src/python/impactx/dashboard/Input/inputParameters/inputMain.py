"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardBase, CardComponents, InputComponents, generalFunctions
from . import InputFunctions

server, state, ctrl = setup_server()
from ..defaults import TRACKING_MODE_PROPERTIES


class InputParameters(CardBase):
    """
    User-Input section for beam properties.
    """

    HEADER_NAME = "Input Parameters"

    def __init__(self):
        super().__init__()

    @state.change("kin_energy_unit")
    def on_kin_energy_unit_change(**kwargs) -> None:
        if state.kin_energy_on_ui != 0:
            InputFunctions.update_kin_energy_sim_value()

    @state.change("tracking_mode")
    def on_tracking_mode_change(**kwargs) -> None:
        """
        Sync the relevant UI components whenever
        the user selects a new tracking mode.
        """
        ui_props = TRACKING_MODE_PROPERTIES[state.tracking_mode]
        for prop_name, prop in ui_props.items():
            setattr(state, prop_name, prop)

        current_sc_list = ui_props.get("space_charge_list", [])
        if state.space_charge not in current_sc_list:
            state.space_charge = current_sc_list[0]
        generalFunctions.update_simulation_validation_status()

    def card_content(self):
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(self.HEADER_NAME)
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=12):
                        InputComponents.select(
                            label="Tracking Mode",
                        )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(md=6, sm=12):
                        InputComponents.select(
                            label="Space Charge",
                            items=("space_charge_list",),
                            disabled=("disable_space_charge",),
                        )
                    with vuetify.VCol(
                        cols=3, classes="py-0 d-flex justify-center align-center"
                    ):
                        InputComponents.checkbox(
                            label="CSR",
                            disabled=("disable_csr",),
                        )
                    with vuetify.VCol(
                        cols=3, classes="py-0 d-flex justify-center align-center"
                    ):
                        InputComponents.checkbox(
                            label="ISR",
                            disabled=("disable_isr",),
                        )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=6):
                        InputComponents.text_field(
                            label="Ref. Particle Charge",
                            v_model_name="charge_qe",
                        )
                    with vuetify.VCol(cols=6):
                        InputComponents.text_field(
                            label="Ref. Particle Mass",
                            v_model_name="mass_MeV",
                        )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=12):
                        InputComponents.text_field(
                            label="Number of Particles",
                            v_model_name="npart",
                        )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=8):
                        InputComponents.text_field(
                            label="Kinetic Energy",
                            v_model_name="kin_energy_on_ui",
                        )
                    with vuetify.VCol(cols=4):
                        InputComponents.select(
                            label="Unit",
                            v_model_name="kin_energy_unit",
                        )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=12):
                        InputComponents.text_field(
                            label="Bunch Charge",
                            v_model_name="bunch_charge_C",
                        )
