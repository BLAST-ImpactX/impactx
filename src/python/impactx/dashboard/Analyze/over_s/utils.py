"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""
import glob
import os
from ... import setup_server
from . import line_plot_1d
from ..analyzeFunctions import AnalyzeFunctions
server, state, ctrl = setup_server()

DEFAULT_HEADERS = ["s", "beta_x", "beta_y"]

state.selected_headers = DEFAULT_HEADERS
state.filtered_data = []
state.all_data = []
state.all_headers = []

@state.change("selected_headers")
def on_header_selection_change(selected_headers, **kwargs):
    state.filtered_headers = AnalyzeFunctions.filter_headers(
        state.all_headers, selected_headers
    )
    state.filtered_data = AnalyzeFunctions.filter_data(state.all_data, selected_headers)

@state.change("filtered_data", "active_plot")
def on_filtered_data_change(**kwargs):
    over_s.update()

class VisualizeOverS:
    def _update_table(self):
        """
        Combines reducedBeam and refParticle files
        and updates data table upon column selection by user
        """

        self.load_dataTable_data()
        state.filtered_data = AnalyzeFunctions.filter_data(
            state.all_data, state.selected_headers
        )
        state.filtered_headers = AnalyzeFunctions.filter_headers(
            state.all_headers, state.selected_headers
        )

    def _update_plot(self):
        """
        Updates the plot in the 'Plot Over S' tab
        """

        fig = line_plot_1d(state.selected_headers, state.filtered_data)
        ctrl.plotly_figure_update(fig)


    def update(self):
        """
        Updates the 'Plot Over S' tab with the latest data and plot.
        Called once when the simulation is complete.
        """
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
        state.headers_without_step_or_s = state.all_headers[2:]

over_s = VisualizeOverS()