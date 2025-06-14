"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import importlib

import pytest
from seleniumbase import SB

from .utils import (
    DashboardTester,
    start_dashboard,
    wait_for_interaction_ready,
    wait_for_server_ready,
)


@pytest.mark.skipif(
    importlib.util.find_spec("seleniumbase") is None,
    reason="seleniumbase is not available",
)
def test_dashboard():

    """
    Tests the ImpactX dashboard with inputs from 'examples/fodo/run_fodo.py'.

    The test includes:
    - Launching the dashboard
    - Configuring the beam properties
    - Configuring the beam distribution
    - Configuring the lattice
        - Through direct numeric input
        - Through the variable configuration
    - Running the simulation
    - Checking if the simulation completes successfully
    """
    app_process = None

    try:
        with SB(headless=True) as sb:
            dashboard = DashboardTester(sb)

            # Setup Dashboard
            app_process = start_dashboard()
            wait_for_server_ready(app_process)
            sb.open("http://localhost:8080/index.html#/Input")
            wait_for_interaction_ready(sb)

            # Adjust beam properties
            BEAM_PROPERTIES = {
                "tracking_mode": "Particle Tracking",
                "space_charge": "false",
                "csr": False,
                "isr": False,
                "charge_qe": -1,
                "mass_MeV": 0.510998950,
                "npart": 10000,
                "bunch_charge_C": 1e-9,
            }

            for param_id, param_value in BEAM_PROPERTIES.items():
                dashboard.set_input(param_id, param_value)

            # Adjust beam distribution
            DISTRIBUTION_PARAMETERS = {
                "distribution": "Waterbag",
                "distribution_type": "Quadratic",
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
                dashboard.set_input(param_id, param_value)

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
                dashboard.set_input(param_id, param_value)

            # Add variable
            sb.click("#lattice_settings")
            VARIABLES = {
                "variable_name_1": "ns",
                "variable_value_1": 25,
            }

            for param_id, param_value in VARIABLES.items():
                dashboard.set_input(param_id, param_value)

            # Run the simulation and check if it completes successfully
            sb.click("#Run_route")
            sb.click("#run_simulation_button")
            dashboard.assert_state("sim_progress_status", "Complete!")

    finally:
        if app_process is not None:
            app_process.terminate()
