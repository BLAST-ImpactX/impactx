"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import glob
import os

import pandas as pd

from ... import ctrl, state
from ...Run.utils import SimulationHelper
from .plot import over_s_plot

DEFAULT_HEADERS = ["beta_x", "beta_y"]
UNSELECTABLE_HEADERS = ["step", "s"]

state.selected_headers = DEFAULT_HEADERS
state.over_s_possible_data = []
state.over_s_possible_headers = []


@state.change("over_s_possible_headers")
def on_over_s_possible_headers_updated(**_):
    state.selectable_headers = [
        header
        for header in state.over_s_possible_headers
        if header["key"] not in UNSELECTABLE_HEADERS
    ]


@state.change("selected_headers")
def on_header_selection_change(**_):
    """
    Called whenever the selected headers
    for the 'Plot Over S' visualization is changed.
    """
    over_s._update_table()


@state.change("over_s_table_headers")
def on_over_s_table_headers_change(**_):
    """
    Called whenever the data table for the 'Plot Over S' visualization is changed.
    """
    over_s._update_plot()


class VisualizeOverS:
    def update(self):
        """
        Updates the 'Plot Over S' tab with the latest data and plot.
        Called once when the simulation is complete.
        """
        self._load_data_table()
        self._update_table()
        self._update_plot()

    def _update_table(self):
        """
        Updates the data table for the 'Plot Over S' visualization.
        """

        # Only display the headers that are selected by the user
        selected_headers = set(state.selected_headers)
        state.over_s_table_headers = [
            header
            for header in state.over_s_possible_headers
            if header["key"] in selected_headers
        ]

        # Display the corresponding data rows
        state.over_s_table_data = [
            row
            for row in state.over_s_possible_data
            if any(key in selected_headers for key in row.keys())
        ]

    def _update_plot(self):
        """
        Updates the plot for the 'Plot Over S' visualization.
        """
        ctrl.plotly_figure_update(over_s_plot())

    def _load_data_table(self) -> None:
        """
        When called, retrieves both beam and reference particle data, combines
        them into a single DataFrame, and updates
        """

        DIAGS_DIRECTORY = os.path.join(os.getcwd(), "diags")
        REDUCED_BEAM_FILE = glob.glob(
            f"{DIAGS_DIRECTORY}/reduced_beam_characteristics.*"
        )
        REF_PARTICLE_FILE = glob.glob(f"{DIAGS_DIRECTORY}/ref_particle.*")

        if not REDUCED_BEAM_FILE or not REF_PARTICLE_FILE:
            return

        beam_data = pd.read_csv(REDUCED_BEAM_FILE[0], sep=" ")
        ref_particle_data = pd.read_csv(REF_PARTICLE_FILE[0], sep=" ")
        combined_data = pd.merge(beam_data, ref_particle_data, how="outer")

        over_s._update_possible_data(combined_data)

    def _update_possible_data(self, data: pd.DataFrame) -> None:
        """
        Updates the possible headers and data for the 'Plot Over S' visualization.
        """
        state.over_s_possible_headers = [
            {"title": header, "key": header} for header in data.columns
        ]
        state.over_s_possible_data = data.to_dict(orient="records")

        SimulationHelper.save_over_s_table_output()


over_s = VisualizeOverS()
