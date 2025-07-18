ERROR_INPUT = "error"
NON_ERROR_INPUT = "no_error"

def test_validation(dashboard):
    """
    Test the errors_tracker and validation of the dashboard.
    We set strings to all Non-Dropdown inputs and verify the count of tracked errors.
    """

    SIMULATION = {
        "tracking_mode": "Particle Tracking",
        "csr": True,
        "space_charge": "3D",
        "charge_qe": ERROR_INPUT,
        "mass_MeV": ERROR_INPUT,
        "kin_energy_on_ui": ERROR_INPUT,
        "npart": ERROR_INPUT,
        "bunch_charge_C": ERROR_INPUT,
    }

    CSR = {
        "csr_bins": ERROR_INPUT,
    }

    SPACE_CHARGE = {
        "max_level": 4,
        "poisson_solver": "multigrid",
        "n_cell_x": ERROR_INPUT,
        "n_cell_y": ERROR_INPUT,
        "n_cell_z": ERROR_INPUT,
        "blocking_factor_x": ERROR_INPUT,
        "blocking_factor_y": ERROR_INPUT,
        "blocking_factor_z": ERROR_INPUT,
        "prob_relative_0": ERROR_INPUT,
        "prob_relative_1": ERROR_INPUT,
        "prob_relative_2": ERROR_INPUT,
        "prob_relative_3": ERROR_INPUT,
        "prob_relative_4": ERROR_INPUT,
        "mlmg_relative_tolerance": ERROR_INPUT,
        "mlmg_absolute_tolerance": ERROR_INPUT,
        "mlmg_verbosity": ERROR_INPUT,
        "mlmg_max_iters": ERROR_INPUT,
    }

    DISTRIBUTION_PARAMETERS = {
        "distribution_type": "Quadratic",
        "lambdaX": ERROR_INPUT,
        "lambdaY": ERROR_INPUT,
        "lambdaT": ERROR_INPUT,
        "lambdaPx": ERROR_INPUT,
        "lambdaPy": ERROR_INPUT,
        "lambdaPt": ERROR_INPUT,
        "muxpx": ERROR_INPUT,
        "muypy": ERROR_INPUT,
        "mutpt": ERROR_INPUT,
    }

    INPUTS = SIMULATION  | DISTRIBUTION_PARAMETERS | SPACE_CHARGE | CSR
    for param_id, value in INPUTS.items():
        dashboard.set_input(param_id, value)


    dashboard.sb.click("#multigrid_settings_button")
    MULTIGRID_ADVANCED_SETTINGS = {
        "mlmg_relative_tolerance": ERROR_INPUT,
        "mlmg_absolute_tolerance": ERROR_INPUT,
        "mlmg_verbosity": ERROR_INPUT,
        "mlmg_max_iters": ERROR_INPUT,
    }
    for param_id, value in MULTIGRID_ADVANCED_SETTINGS.items():
        dashboard.set_input(param_id, value)

    # Lattice
    dashboard.add_lattice_element("CFbend")
    LATTICE_PARAMS = {
        "ds1": ERROR_INPUT,
        "rc1": ERROR_INPUT,
        "k1": ERROR_INPUT,
        "name1": NON_ERROR_INPUT,
    }

    for param_id, value in LATTICE_PARAMS.items():
        dashboard.set_input(param_id, value)

    # Variable
    dashboard.sb.click("#lattice_settings")
    VARIABLES = {"variable_name_1": NON_ERROR_INPUT, "variable_value_1": ERROR_INPUT}
    for param_id, value in VARIABLES.items():
        dashboard.set_input(param_id, value)

    # Check to make sure that the input_errors is only tracking what is shown on the UI
    assert dashboard.get_state("number_of_input_errors") == 33
    dashboard.set_input("poisson_solver", "fft")
    assert dashboard.get_state("number_of_input_errors") == 29
    dashboard.set_input("csr", False)
    dashboard.set_input("space_charge", "false")
    assert dashboard.get_state("number_of_input_errors") == 18
    dashboard.set_input("reset_lattice_button", True)
    assert dashboard.get_state("number_of_input_errors") == 15