"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import setup_server, vuetify

server, state, ctrl = setup_server()


class LatticeVisualizerDialogs:
    @staticmethod
    def element_colors_tab():
        with vuetify.VCardText():
            with vuetify.VRow():
                with vuetify.VCol(cols=12):
                    vuetify.VCardSubtitle(
                        "Element Color Mapping",
                    )

    @staticmethod
    def general_settings_tab():
        with vuetify.VCardText():
            with vuetify.VRow():
                with vuetify.VCol(cols=12):
                    vuetify.VCardSubtitle(
                        "Settings",
                    )
