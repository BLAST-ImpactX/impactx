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

    # Wait for the file to finish loading and processing
    # The importing_file state is True while loading, False when done
    dashboard.assert_state("importing_file", False)

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
        ("#ds1", 0.25),
        ("#ds2", 1),
        ("#k2", 1.0),
        ("#ds3", 0.5),
        ("#ds4", 1.0),
        ("#k4", -1.0),
        ("#ds5", 0.25),
    ]

    # Check state parameters
    for param_name, expected_value in (
        BEAM_PARAMETERS + DISTRIBUTION_PARAMETERS + LATTICE_INPUTS
    ):
        dashboard.assert_state(param_name, expected_value)

    # Check input values
    for element_id, expected_value in DISTRIBUTION_VALUES + LATTICE_CONFIGURATION:
        # Wait for element to be present in the DOM (doesn't require visibility)
        # This is important for CI environments where rendering may be slower
        dashboard.sb.wait_for_element_present(element_id, timeout=10)

        # Wait for the element's value to be populated
        # DOM elements may take time to reflect updated state values after file loading
        dashboard.wait_for_element_value(element_id, timeout=10)

        # Get the value attribute - this works even if the element isn't visible
        actual_value = float(dashboard.sb.get_attribute(element_id, "value"))
        assert actual_value == pytest.approx(expected_value, **APPROX_TOL), (
            f"{element_id}: expected {expected_value}, got {actual_value}"
        )
