"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from ... import setup_server

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

CONVERSION_FACTORS = {
    "meV": 1.0e-9,
    "eV": 1.0e-6,
    "keV": 1.0e-3,
    "MeV": 1.0,
    "GeV": 1.0e3,
    "TeV": 1.0e6,
}


class InputFunctions:
    """
    Helper functions for the beam properties.
    """

    @staticmethod
    def update_kin_energy_sim_value() -> None:
        """
        Converts UI input to MeV for internal simulation value.
        """

        state.kin_energy_MeV = (
            state.kin_energy_on_ui * CONVERSION_FACTORS[state.kin_energy_unit]
        )
