def test_dashboard(dashboard):
    """
    End-to-end test of the ImpactX dashboard by directly setting inputs via UI,
    instead of importing a file. Mirrors the configuration in 'examples/fodo/run_fodo.py'.

    The test includes:
    - Configuring the beam properties / beam distribution
    - Configuring the lattice
        - Through direct numeric input
        - Through the variable configuration
    - Checking if the lattice statistics update correctly
    - Running the simulation
    - Checking if the simulation completes successfully
    """

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

    for param_id, value in {**BEAM_PROPERTIES, **DISTRIBUTION_PARAMETERS}.items():
        dashboard.set_input(param_id, value)

    # Lattice
    for name in ["Drift", "Quad", "Drift", "Quad", "Drift"]:
        dashboard.add_lattice_element(name)

    LATTICE_PARAMS = {
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
    for param_id, value in LATTICE_PARAMS.items():
        dashboard.set_input(param_id, value)

    # Variable
    dashboard.sb.click("#lattice_settings")
    VARIABLES = {"variable_name_1": "ns", "variable_value_1": 25}
    for param_id, value in VARIABLES.items():
        dashboard.set_input(param_id, value)

    # Check if the lattice statistics update correctly
    EXPECTED_LATTICE_STATS = {
        "total_elements": 5,
        "total_length": "3.0m",
        "max_length": "1.0m",
        "avg_length": "0.6m",
        "min_length": "0.25m",
        "total_steps": 125,
        "periods": 1,
        "element_counts": {"drift": 3, "quad": 2},
    }

    for state_name, expected_value in EXPECTED_LATTICE_STATS.items():
        dashboard.assert_state(state_name, expected_value)

    # Run simulation
    dashboard.sb.click("#Run_route")
    dashboard.sb.click("#run_simulation_button")
    dashboard.assert_state("sim_progress_status", "Complete!")
