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


TRACKING_MODE_UI_PROPERTIES: dict[str, dict[str, bool]] = {
    "Reference Tracking": {
        "space_charge":        False,
        "csr":                 False,
        "disable_space_charge": True,
        "disable_csr":          True,
    },
    "Envelope Tracking": {
        "csr":                 False,
        "disable_space_charge": False,
        "disable_csr":          True,
    },
    "Particle Tracking": {
        "disable_space_charge": False,
        "disable_csr":          False,
    },
}

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
        ui_props = TRACKING_MODE_UI_PROPERTIES[state.tracking_mode]
        for prop_name, prop in ui_props.items():
            setattr(state, prop_name, prop)

    def card_content(self):
        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(self.HEADER_NAME)
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(md=6, sm=12):
                        InputComponents.select(
                            label="Tracking Mode",
                        )
                    with vuetify.VCol(classes="ga-4 py-0 d-flex align-center"):
                        InputComponents.checkbox(
                            label="SC",
                            v_model_name="space_charge",
                            disabled=("disable_space_charge",),
                        )
                        InputComponents.checkbox(
                            label="CSR",
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
