"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import sys

from trame.ui.router import RouterViewLayout
from trame.ui.vuetify import SinglePageWithDrawerLayout
from trame.widgets import router, xterm

from . import (
    AnalyzeSimulation,
    DistributionParameters,
    GeneralToolbar,
    InputParameters,
    LatticeConfiguration,
    NavigationComponents,
    SpaceChargeConfiguration,
    csrConfiguration,
    setup_server,
    vuetify,
)
from .start import main

server, state, ctrl = setup_server()

from .Input.shared import SharedUtilities

shared_utilities = SharedUtilities()

inputParameters = InputParameters()

with RouterViewLayout(server, "/Input"):
    with vuetify.VContainer(fluid=True):
        with vuetify.VRow():
            with vuetify.VCol(cols="auto", classes="pa-2"):
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        inputParameters.card()
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        SpaceChargeConfiguration.card()
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        csrConfiguration.card()
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        DistributionParameters.card()
                with vuetify.VRow(no_gutters=True):
                    with vuetify.VCol(cols="auto", classes="pa-2"):
                        LatticeConfiguration.card()

with RouterViewLayout(server, "/Analyze"):
    with vuetify.VContainer(fluid=True):
        with vuetify.VRow(no_gutters=True, classes="fill-height"):
            with vuetify.VCol(cols="auto", classes="pa-2 fill-height"):
                AnalyzeSimulation.card()
            with vuetify.VCol(cols="auto", classes="pa-2 fill-height"):
                AnalyzeSimulation.plot()


# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------
def init_terminal():
    with xterm.XTerm(v_if="$route.path == '/Run'") as term:
        ctrl.terminal_print = term.writeln


def application():
    init_terminal()
    with SinglePageWithDrawerLayout(server) as layout:
        layout.title.hide()
        with layout.toolbar:
            with vuetify.Template(v_if="$route.path == '/Analyze'"):
                GeneralToolbar.dashboard_toolbar("analyze")
            with vuetify.Template(v_if="$route.path == '/Input'"):
                GeneralToolbar.dashboard_toolbar("input")
            with vuetify.Template(v_if="$route.path == '/Run'"):
                GeneralToolbar.dashboard_toolbar("run")

        with layout.drawer as drawer:
            drawer.width = 200
            with vuetify.VList():
                vuetify.VSubheader("Simulation")
            NavigationComponents.create_route("Input", "mdi-file-edit")
            NavigationComponents.create_route("Run", "mdi-play")
            NavigationComponents.create_route("Analyze", "mdi-chart-box-multiple")

        with layout.content:
            router.RouterView()
            init_terminal()
    return layout


application()
# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit(main())
