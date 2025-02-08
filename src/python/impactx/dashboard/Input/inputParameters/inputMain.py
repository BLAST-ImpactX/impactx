"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from .. import CardComponents, InputComponents
from . import InputFunctions

server, state, ctrl = setup_server()


class InputParameters:
    """
    User-Input section for beam properties.
    """

    @state.change("kin_energy_unit")
    def on_kin_energy_unit_change(**kwargs) -> None:
        if state.kin_energy_on_ui != 0:
            InputFunctions.update_kin_energy_sim_value()

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
