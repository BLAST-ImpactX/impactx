"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from impactx.impactx_pybind import ImpactX, RefPart

from .. import setup_server
from .defaults_helper import InputDefaultsHelper

server, state, ctrl = setup_server()


class DashboardDefaults:
    """
    Defaults for input parameters in the ImpactX dashboard.
    """

    COLLAPSABLE_SECTIONS = [
        "collapse_input_parameters",
        "collapse_csr",
        "collapse_distribution_parameters",
        "collapse_space_charge",
        "collapse_lattice_configuration",
    ]
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
        "kin_energy_on_ui": 2e3,
        "kin_energy_MeV": 2e3,
        "kin_energy_unit": "MeV",
        "bunch_charge_C": 1e-9,
    }

    DISTRIBUTION_PARAMETERS = {
        "distribution": "Waterbag",
        "distribution_type": "Twiss",
    }

    LATTICE_CONFIGURATION = {
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
        **DISTRIBUTION_PARAMETERS,
        **LATTICE_CONFIGURATION,
        **SPACE_CHARGE,
        **CSR,
        **LISTS,
    }

    TYPES = {
        "npart": "int",
        "kin_energy_on_ui": "float",
        "bunch_charge_C": "float",
        "mass_MeV": "float",
        "charge_qe": "int",
        "csr_bins": "int",
    }

    VALIDATION_CONDITION = {
        "charge_qe": ["non_zero"],
        "mass_MeV": ["positive"],
        "csr_bins": ["positive"],
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

    DOCUMENTATION = {
        "input_parameters": "https://impactx.readthedocs.io/en/latest/usage/python.html#impactx.ImpactX",
        "lattice_configuration": "https://impactx.readthedocs.io/en/latest/usage/python.html#lattice-elements",
        "distribution_parameters": "https://impactx.readthedocs.io/en/latest/usage/python.html#initial-beam-distributions",
        "space_charge": "https://impactx.readthedocs.io/en/latest/usage/parameters.html#space-charge",
        "csr": "https://impactx.readthedocs.io/en/latest/usage/parameters.html#coherent-synchrotron-radiation-csr",
    }


class TooltipDefaults:
    """
    Defaults for input toolips in the ImpactX dashboard.
    """

    TOOLTIP = InputDefaultsHelper.get_docstrings(
        [RefPart, ImpactX], DashboardDefaults.DEFAULT_VALUES
    )


class ToolbarDefaults:
    """
    Default styling and states for the toolbar
    section in the ImpactX dashboard.
    """

    TOOLBAR_SIZE = 64
    FOOTER_SIZE = 8


class UIDefaults:
    """
    Default UI which the input cards reply on in the ImpactX dashboard.
    """

    ROW_STYLE = {
        "dense": False,
    }

    CARD_TEXT_OVERFLOW = {
        "classes": "custom-scrollbar",
        "style": {
            "flex": "1",
            "overflow-y": "auto",
            "overflow-x": "auto",
        },
    }

    CARD_STYLE = {
        "display": "flex",
        "flex-direction": "column",
    }
