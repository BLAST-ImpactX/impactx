"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import CardComponents
from .... import setup_server, vuetify

server, state, ctrl = setup_server()


class LatticeVisualizerDialogs:

    @staticmethod
    def settings():
        """
        A button to open the settings dialog for the lattice visualizer.
        """
        CardComponents.card_button(
            "mdi-cog",
            color="white",
            click="lattice_visualizer_dialog_settings = true",
            description="Settings",
        )

    @staticmethod
    def element_colors_tab():
        """
        A tab inside of the settings dialog to manage element colors.
        """

        with vuetify.VCardText():
            with vuetify.VRow():
                with vuetify.VCol(cols=12):
                    vuetify.VCardSubtitle(
                        "Element Color Mapping",
                    )

    @staticmethod
    def general_settings_tab():
        """
        A tab inside of the settings dialog for general settings.
        """

        with vuetify.VCardText():
            with vuetify.VRow():
                with vuetify.VCol(cols=12):
                    vuetify.VCardSubtitle(
                        "Settings",
                    )
