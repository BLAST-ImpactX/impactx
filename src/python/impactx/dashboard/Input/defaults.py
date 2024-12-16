class DashboardDefaults:
    """
    Defaults for input parameters in the ImpactX dashboard.
    """

    PARAMETER_VALUES = {
        # Input Parameters
        "charge_qe": -1,
        "mass_MeV": 0.51099895,
        "npart": 1000,
        "kin_energy": 2e3,
        "bunch_charge_C": 1e-9,
        "particle_shape": 2,
        # Space Charge
        "poisson_solver": "fft",
        "max_level": 0,
        "n_cell": 32,
        "blocking_factor": 16,
        "prob_relative_first_value_fft": 1.1,
        "prob_relative_first_value_multigrid": 3.1,
        "mlmg_relative_tolerance": 1.0e-7,
        "mlmg_absolute_tolerance": 0,
        "mlmg_verbosity": 1,
        "mlmg_max_iters": 100,
        # CSR
        "csr_bins": 150,
    }

    # If parameter is not included in the dictionary, default step amount is 1.
    PARAMETER_STEPS = {
        # Single input
        "mass_MeV": 0.1,
        "bunch_charge_C": 1e-11,
        "prob_relative": 0.1,
        "mlmg_relative_tolerance": 1e-12,
        "mlmg_absolute_tolerance": 1e-12,
        # Shared inputs (x,y,z)
        "beta": 0.1,
        "emitt": 1e-7,
        "alpha": 0.1,
    }

    PARAMETER_UNITS = {
        # Single input
        "charge_qe": "qe",
        "mass_MeV": "MeV",
        "bunch_charge_C": "C",
        "mlmg_absolute_tolerance": "V/m",
        # Shared inputs (x,y,z)
        "beta": "m",
        "emitt": "m",
    }
