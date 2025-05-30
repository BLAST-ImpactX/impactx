"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import CardComponents


class LatticeVisualizerComponents:
    @staticmethod
    def settings():
        CardComponents.card_button(
            "mdi-cog",
            color="grey-darken-2",
            click="lattice_visualizer_dialog_settings = true",
            documentation="Settings",
        )
