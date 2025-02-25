"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardBase, CardComponents, InputComponents
from . import InputFunctions

server, state, ctrl = setup_server()


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
        if state.tracking_mode in ("Envelope Tracking", "Reference Tracking"):
            state.space_charge = False
            state.csr = False

            state.disable_space_charge = True
            state.disable_csr = True
        else:
            state.disable_space_charge = False
            state.disable_csr = False

    def card_content(self):
        with vuetify.VCard(style=self.collapsable):
            CardComponents.input_header(self.HEADER_NAME)
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=5.5):
                        InputComponents.select(
                            label="Tracking Mode",
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VCheckbox(
                            label="Space Charge",
                            v_model=("space_charge", False),
                            dense=True,
                            disabled=("disable_space_charge",),
                        )
                    with vuetify.VCol(cols="auto"):
                        vuetify.VCheckbox(
                            label="CSR",
                            v_model=("csr", False),
                            dense=True,
                            disabled=("disable_csr",),
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
