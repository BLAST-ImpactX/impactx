class DashboardDefaults:
    """
    Defaults for input parameters in the ImpactX dashboard.
    """

    # -------------------------------------------------------------------------
    # Inputs by section
    # -------------------------------------------------------------------------

    SELECTION = {
        "space_charge": False,
        "csr": False,
    }

    INPUT_PARAMETERS = {
        "charge_qe": -1,
        "mass_MeV": 0.51099895,
        "npart": 1000,
        "kin_energy": 2e3,
        "kin_energy_unit": "MeV",
        "bunch_charge_C": 1e-9,
    }

    DISTRIBUTION = {
        "selected_distribution": "Waterbag",
        "selected_distribution_type": "Twiss",
    }

    LATTICE = {
        "selected_lattice_list": [],
        "selected_lattice": None,
    }

    SPACE_CHARGE = {
        "dynamic_size": False,
        "poisson_solver": "fft",
        "particle_shape": 2,
        "max_level": 0,
        "n_cell_x": 32,
        "n_cell_y": 32,
        "n_cell_z": 32,
        "blocking_factor_x": 16,
        "blocking_factor_y": 16,
        "blocking_factor_z": 16,
        "prob_relative_first_value_fft": 1.1,
        "prob_relative_first_value_multigrid": 3.1,
        "mlmg_relative_tolerance": 1.0e-7,
        "mlmg_absolute_tolerance": 0,
        "mlmg_verbosity": 1,
        "mlmg_max_iters": 100,
    }

    CSR = {
        "particle_shape": 2,
        "csr_bins": 150,
    }

    LISTS = {
        "kin_energy_unit_list": ["meV", "eV", "keV", "MeV", "GeV", "TeV"],
        "distribution_type_list": ["Twiss", "Quadratic"],
        "poisson_solver_list": ["fft", "multigrid"],
        "particle_shape_list": [1, 2, 3],
        "max_level_list": [0, 1, 2, 3, 4],
    }

    # -------------------------------------------------------------------------
    # Main
    # -------------------------------------------------------------------------

    DEFAULT_VALUES = {
        **SELECTION,
        **INPUT_PARAMETERS,
        **DISTRIBUTION,
        **LATTICE,
        **SPACE_CHARGE,
        **CSR,
        **LISTS,
    }

    # If a parameter is not included in the dictionary, default step amount is 1.
    STEPS = {
        "mass_MeV": 0.1,
        "bunch_charge_C": 1e-11,
        "prob_relative": 0.1,
        "mlmg_relative_tolerance": 1e-12,
        "mlmg_absolute_tolerance": 1e-12,
        "beta": 0.1,
        "emitt": 1e-7,
        "alpha": 0.1,
    }

    UNITS = {
        "charge_qe": "qe",
        "mass_MeV": "MeV",
        "bunch_charge_C": "C",
        "mlmg_absolute_tolerance": "V/m",
        "beta": "m",
        "emitt": "m",
    }
