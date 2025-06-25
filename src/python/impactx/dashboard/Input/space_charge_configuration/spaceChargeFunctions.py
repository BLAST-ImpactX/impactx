"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server
from .. import CardComponents

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


class SpaceChargeFunctions:
    @staticmethod
    def multigrid_settings():
        CardComponents.card_button(
            "mdi-cog",
            color="grey-darken-2",
            click="space_charge_dialog_settings = true",
            v_if="poisson_solver == 'multigrid'",
            documentation="Settings",
        )
