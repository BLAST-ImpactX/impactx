"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""
import glob
import os

from . import line_plot_1d
from ..analyzeFunctions import AnalyzeFunctions
from ... import setup_server

server, state, ctrl = setup_server()

DEFAULT_HEADERS = ["s", "beta_x", "beta_y"]

state.selected_headers = DEFAULT_HEADERS
state.all_data = []
state.all_headers = []

@state.change("selected_headers")
def on_header_selection_change(**_):
    over_s.update()

class VisualizeOverS:
    def _update_table(self):
        """
        Combines reducedBeam and refParticle files
        and updates data table upon column selection by user
        """

        state.over_s_data = AnalyzeFunctions.filter_data()
        state.over_s_headers = AnalyzeFunctions.filter_headers()

    def _update_plot(self):
        """
        Updates the plot in the 'Plot Over S' tab
        """

        fig = line_plot_1d(state.selected_headers, state.over_s_data)
        ctrl.plotly_figure_update(fig)


    def update(self):
        """
        Updates the 'Plot Over S' tab with the latest data and plot.
        Called once when the simulation is complete.
        """
        self.load_dataTable_data()
        self._update_table()
        self._update_plot()

    def load_dataTable_data(self):
        """
        Loads and processes data from combined beam and reference particle files.
        """

        CURRENT_DIR = os.getcwd()
        DIAGS_DIR = os.path.join(CURRENT_DIR, "diags")

        base_path = DIAGS_DIR + "/"
        REDUCED_BEAM_DATA = glob.glob(base_path + "reduced_beam_characteristics.*")[0]
        REF_PARTICLE_DATA = glob.glob(base_path + "ref_particle.*")[0]

        if not os.path.exists(REDUCED_BEAM_DATA) or not os.path.exists(REF_PARTICLE_DATA):
            ctrl.terminal_print(
                "Diagnostics files are missing. Please ensure they are in the correct directory."
            )
            return

        combined_files = AnalyzeFunctions.combine_files(
            REDUCED_BEAM_DATA, REF_PARTICLE_DATA
        )
        combined_files_data_converted_to_dictionary_format = (
            AnalyzeFunctions.convert_to_dict(combined_files)
        )
        data, headers = combined_files_data_converted_to_dictionary_format
        state.all_data = data
        state.all_headers = headers
        state.selectable_headers = [
            header for header in state.all_headers if header["key"] not in ("step", "s")
        ]

over_s = VisualizeOverS()