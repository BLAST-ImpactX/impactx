"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import time

import pytest

from impactx import RefPart

from .utils import APPROX_TOL, TIMEOUT


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
        ("mass_MeV", RefPart().set_species("electron").mass_MeV),
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

    # Check input values using get_attribute (does not require element
    # visibility, unlike get_value which fails for scrollable containers).
    for element_id, expected_value in DISTRIBUTION_VALUES + LATTICE_CONFIGURATION:
        for i in range(TIMEOUT):
            actual_value = dashboard.sb.get_attribute(element_id, "value")
            if actual_value is not None and actual_value != "":
                break
            time.sleep(1)
        else:
            raise AssertionError(
                f"{element_id}: value was still empty after {TIMEOUT} seconds"
            )

        if isinstance(expected_value, str):
            assert actual_value == expected_value, (
                f"{element_id}: expected '{expected_value}', got '{actual_value}'"
            )
        else:
            actual_value = float(actual_value)
            assert actual_value == pytest.approx(expected_value, **APPROX_TOL), (
                f"{element_id}: expected {expected_value}, got {actual_value}"
            )


def test_python_import_envelope(dashboard):
    """
    Imports a proton envelope with Twiss parameters and a repeated lattice cell
    from 'testdata/envelope_proton.py' and runs it.
    """
    dashboard.load_example("testdata/envelope_proton.py", manual=True)

    BEAM_PARAMETERS = [
        ("tracking_mode", "Envelope Tracking"),
        ("space_charge", "3D"),
        ("charge_qe", 1),
        ("mass_MeV", RefPart().set_species("proton").mass_MeV),
        ("bunch_charge_C", 1e-12),
        ("distribution_type", "Twiss"),
        ("total_elements", 9),
        ("total_length", "1.3m"),
    ]
    for param_name, expected_value in BEAM_PARAMETERS:
        dashboard.assert_state(param_name, expected_value)

    lattice = dashboard.get_state("selected_lattice_list")
    for index, name, value in ((0, "name", "start"), (1, "modify_ref_part", False)):
        parameter = next(
            p for p in lattice[index]["parameters"] if p["parameter_name"] == name
        )
        assert parameter["sim_input"] == value, parameter

    TWISS_PARAMETERS = {"emitt_x": 1e-6, "emitt_y": 2e-6, "alpha_x": -0.5}
    for i in range(TIMEOUT):
        parameters = dashboard.get_state("selected_distribution_parameters") or {}
        values = {name: parameters.get(name, {}).get("value") for name in parameters}
        if all(
            isinstance(values.get(name), (int, float))
            and values[name] == pytest.approx(expected, **APPROX_TOL)
            for name, expected in TWISS_PARAMETERS.items()
        ):
            break
        time.sleep(1)
    else:
        raise AssertionError(f"Twiss parameters {values} != {TWISS_PARAMETERS}")

    dashboard.sb.click("#Run_route")
    dashboard.sb.click("#run_simulation_button")
    dashboard.assert_state("sim_progress_status", "Complete!")
    assert not dashboard.get_state("import_file_error"), dashboard.get_state(
        "import_file_error_message"
    )
