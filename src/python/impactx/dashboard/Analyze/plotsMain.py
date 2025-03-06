"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import base64
import glob
import io
import os

from trame.widgets import matplotlib, plotly

from .. import setup_server, vuetify
from . import AnalyzeFunctions, adjusted_settings_plot, line_plot_1d

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------


def fig_to_base64(fig):
    """
    Puts png in trame-compatible form
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


# Call plot_over_s
def plot_over_s():
    """
    Generates a plot.
    """

    fig = line_plot_1d(state.selected_headers, state.filtered_data)
    ctrl.plotly_figure_update(fig)


def generate_phase_space(pc):
    fig = adjusted_settings_plot(pc)
    ctrl.matplotlib_figure_update(fig)

    fig_original = pc.plot_phasespace()

    if fig_original is not None:
        image_base64 = fig_to_base64(fig_original)
        state.phase_space_png = f"data:image/png;base64, {image_base64}"


PLOTS = {
    "Plot Over S": plot_over_s,
    "Phase Space Plots": None,
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def available_plot_options(simulationClicked):
    """
    Displays plot_options for users based on status of simulation.
    :param simulationClicked (bool): status of simulation status
    :return: list of plot_options for users
    """

    if simulationClicked:
        return list(PLOTS.keys())
    else:
        return ["Run Simulation To See Options"]


def load_dataTable_data():
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


# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------

DEFAULT_HEADERS = ["s", "beta_x", "beta_y"]

state.selected_headers = DEFAULT_HEADERS
state.plot_options = available_plot_options(simulationClicked=False)
state.show_table = False
state.active_plot = None
state.filtered_data = []
state.all_data = []
state.all_headers = []
state.phase_space_png = None

# -----------------------------------------------------------------------------
# Functions to update table/plot
# -----------------------------------------------------------------------------


def update_data_table():
    """
    Combines reducedBeam and refParticle files
    and updates data table upon column selection by user
    """

    load_dataTable_data()
    state.filtered_data = AnalyzeFunctions.filter_data(
        state.all_data, state.selected_headers
    )
    state.filtered_headers = AnalyzeFunctions.filter_headers(
        state.all_headers, state.selected_headers
    )


def update_plot():
    """
    Performs actions to display correct information,
    based on the plot optin selected by the user
    """

    if state.active_plot == "Plot Over S":
        update_data_table()
        plot_over_s()
        state.show_table = True
    elif state.active_plot == "Phase Space Plots":
        state.show_table = False


# -----------------------------------------------------------------------------
# State changes
# -----------------------------------------------------------------------------


@state.change("selected_headers")
def on_header_selection_change(selected_headers, **kwargs):
    state.filtered_headers = AnalyzeFunctions.filter_headers(
        state.all_headers, selected_headers
    )
    state.filtered_data = AnalyzeFunctions.filter_data(state.all_data, selected_headers)


@state.change("filtered_data", "active_plot")
def on_filtered_data_change(**kwargs):
    update_plot()


# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------


class AnalyzeSimulation:
    """
    Prepares contents for the 'Analyze' page.
    """

    @staticmethod
    def card():
        """
        Displays any non-plot content for 'Analyze' page.
        """

        with vuetify.VContainer():
            with vuetify.VCard(v_if=("show_table")):
                with vuetify.VCol(style="width: 500px;"):
                    vuetify.VSelect(
                        v_model=("selected_headers",),
                        items=("headers_without_step_or_s",),
                        label="Select data to view",
                        multiple=True,
                        item_title="title",
                        item_value="key",
                        density="compact",
                        variant="underlined",
                    )
                    vuetify.VDivider()
                    vuetify.VDataTable(
                        headers=("filtered_headers",),
                        items=("filtered_data",),
                        header_class="centered-header",
                        density="compact",
                        height="325px",
                    )

    @staticmethod
    def plot():
        """
        Displays any plot content for 'Analyze' page.
        """

        with vuetify.VContainer(fluid=True):
            with vuetify.VContainer(
                v_if="active_plot === 'Phase Space Plots'", fluid=True
            ):
                with vuetify.VCard(style="height: 50vh; width: 150vh;"):
                    with vuetify.VTabs(v_model=("active_tab", 0)):
                        vuetify.VTab("Plot")
                        vuetify.VTab("Interact")
                    vuetify.VDivider()
                    with vuetify.VTabsWindow(v_model="active_tab"):
                        with vuetify.VTabsWindowItem():
                            vuetify.VImg(
                                v_if=("phase_space_png",), src=("phase_space_png",)
                            )

                        with vuetify.VTabsWindowItem():
                            with vuetify.VContainer(
                                style="height: 37vh; width: 147vh;"
                            ):
                                matplotlib_figure = matplotlib.Figure(
                                    style="position: absolute"
                                )
                                ctrl.matplotlib_figure_update = matplotlib_figure.update

        with vuetify.VContainer(
            v_if="active_plot === 'Plot Over S'", style="height: 50vh; width: 90vh;"
        ):
            plotly_figure = plotly.Figure(
                display_mode_bar="true", config={"responsive": True}
            )
            ctrl.plotly_figure_update = plotly_figure.update
