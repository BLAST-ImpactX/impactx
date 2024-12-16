class DashboardDefaults:
    """
    Defaults for input parameters in the ImpactX dashboard.
    """

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
