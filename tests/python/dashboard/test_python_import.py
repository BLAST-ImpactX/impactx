"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import pytest

from .utils import APPROX_TOL


def test_python_import(dashboard):
    """
    End-to-end test of the ImpactX dashboard by importing input values from a Python file.

    This test mirrors the configuration in 'testdata/example.py' and verifies:
    - Trame state reflects the correct beam and distribution values
    - Form fields (distribution and lattice) are populated with expected inputs

    The configuration is loaded automatically from the file, unlike test_dashboard,
    which sets values manually through direct UI interaction.
    """
    dashboard.load_example("testdata/example.py", manual=True)

    BEAM_PARAMETERS = [
        ("tracking_mode", "Particle Tracking"),
        ("space_charge", "false"),
        ("charge_qe", -1),
        ("mass_MeV", 0.510998950),
        ("npart", 10000),
        ("bunch_charge_C", 1e-9),
    ]

    DISTRIBUTION_PARAMETERS = [
        ("distribution", "Waterbag"),
        ("distribution_type", "Quadratic"),
    ]

    LATTICE_INPUTS = [
        ("periods", 2),
    ]

    DISTRIBUTION_VALUES = [
        ("#lambdaX", 3.9984884770e-5),
        ("#lambdaY", 3.9984884770e-5),
        ("#lambdaT", 1.0e-3),
        ("#lambdaPx", 2.6623538760e-5),
        ("#lambdaPy", 2.6623538760e-5),
        ("#lambdaPt", 2.0e-3),
        ("#muxpx", -0.846574929020762),
        ("#muypy", 0.846574929020762),
        ("#mutpt", 0.0),
    ]

    LATTICE_CONFIGURATION = [
        ("#name1", "drift1"),
        ("#name2", "quad1"),
        ("#name3", "drift2"),
        ("#name4", "quad2"),
        ("#name5", "drift3"),
        ("#name6", "drift3"),
        ("#name7", "quad2"),
        ("#name8", "drift2"),
        ("#ds8", 0.5),
    ]

    # Check state parameters
    for param_name, expected_value in (
        BEAM_PARAMETERS + DISTRIBUTION_PARAMETERS + LATTICE_INPUTS
    ):
        dashboard.assert_state(param_name, expected_value)

    # Check input values
    for element_id, expected_value in DISTRIBUTION_VALUES + LATTICE_CONFIGURATION:
        actual_value = dashboard.sb.get_value(element_id)

        if isinstance(expected_value, str):
            # Compare strings directly
            assert actual_value == expected_value, (
                f"{element_id}: expected '{expected_value}', got '{actual_value}'"
            )
        else:
            # Compare numbers with tolerance
            actual_value = float(actual_value)
            assert actual_value == pytest.approx(expected_value, **APPROX_TOL), (
                f"{element_id}: expected {expected_value}, got {actual_value}"
            )
