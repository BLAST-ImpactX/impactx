"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import importlib

import pytest
from seleniumbase import SB
from utils import (
    DashboardTester,
    start_dashboard,
    wait_for_dashboard,
    wait_for_ready,
)

TIMEOUT = 20


@pytest.mark.skipif(
    importlib.util.find_spec("seleniumbase") is None,
    reason="seleniumbase is not available",
)
def test_dashbnoard():
    """
    ImpactX dashbord testing with inputs from 'examples/fodo/run_fodo.py'.

    Distribution, lattice, and variable inputs are individual states, rather part of a nested state
    which is why we modify those values by the DOM.
    """
    app_process = None

    try:
        with SB(headless=True) as sb:
            dashboard = DashboardTester(sb)

            # Setup Dashboard
            app_process = start_dashboard()
            wait_for_dashboard(app_process, timeout=TIMEOUT)
            sb.open("http://localhost:8080/index.html#/Input")
            wait_for_ready(sb, TIMEOUT)

            # Adjust beam properties
            dashboard.set_state("tracking_mode", "Particle Tracking")
            dashboard.set_state("space_charge", "false")
            # dashboard.set_state("csr", "false")
            # dashboard.set_state("isr", "false")

            dashboard.set_state("charge_qe", -1)
            dashboard.set_state("mass_MeV", 0.510998950)
            dashboard.set_state("npart", 10000)
            dashboard.set_state("bunch_charge_C", 1e-9)

            # Adjust beam distribution
            dashboard.set_state("distribution", "Waterbag")
            dashboard.set_state("distribution_type", "Quadratic")
            DISTRIBUTION_PARAMETERS = {
                "lambdaX": 3.9984884770e-5,
                "lambdaY": 3.9984884770e-5,
                "lambdaT": 1.0e-3,
                "lambdaPx": 2.6623538760e-5,
                "lambdaPy": 2.6623538760e-5,
                "lambdaPt": 2.0e-3,
                "muxpx": -0.846574929020762,
                "muypy": 0.846574929020762,
                "mutpt": 0.0,
            }

            for param_id, param_value in DISTRIBUTION_PARAMETERS.items():
                dashboard.set_js_input(param_id, param_value)

            # Build lattice
            for element_name in ["Drift", "Quad", "Drift", "Quad", "Drift"]:
                dashboard.add_lattice_element(element_name)

            LATTICE_CONFIGURATION = {
                "ds1": 0.25,
                "nslice1": "ns",
                "ds2": 1.0,
                "k2": 1.0,
                "nslice2": "ns",
                "ds3": 0.5,
                "nslice3": "ns",
                "ds4": 1.0,
                "k4": -1.0,
                "nslice4": "ns",
                "ds5": 0.25,
                "nslice5": "ns",
            }

            for param_id, param_value in LATTICE_CONFIGURATION.items():
                dashboard.set_js_input(param_id, param_value)

            # Add variable
            sb.click("#lattice_settings")
            VARIABLES = {
                "variable_name_1": "ns",
                "variable_value_1": 25,
            }

            for param_id, param_value in VARIABLES.items():
                dashboard.set_js_input(param_id, param_value)

            # Run the simulation and check if it completes successfully
            sb.click("#Run_route")
            sb.click("#run_simulation_button")
            dashboard.assert_state("sim_progress_status", "Complete!")

    finally:
        if app_process is not None:
            app_process.terminate()
