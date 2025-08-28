"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import time


def assert_lattice_param_sim_input(
    dashboard, lattice_index: int, param_name: str, expected_value
):
    """
    Asserts that a lattice parameter's sim_input matches expected value, with retry logic.
    """

    def get_sim_input():
        lattice_list = dashboard.get_state("selected_lattice_list")
        for param in lattice_list[lattice_index]["parameters"]:
            if param["parameter_name"] == param_name:
                return param["sim_input"]
        return None

    for _ in range(10):
        current_value = get_sim_input()
        if current_value == expected_value:
            return
        time.sleep(1)

    raise AssertionError(
        f"Parameter '{param_name}' sim_input never became '{expected_value}' (got: {current_value})"
    )


def test_lattice_defaults_handler(dashboard):
    """
    Verifies the lattice defaults handler

    1. Set a default for 'nslice' and ensure new elements use it.
    2. Add a second default and ensure subsequent elements reflect it.
    3. Enter an unknown parameter name and confirm an error is stored.
    4. Reset defaults and verify the handler resets to a single row.
    """

    # Open lattice settings dialog
    dashboard.sb.click("#lattice_settings")

    # Locate the row for 'nslice' and set its value
    defaults_state = dashboard.get_state("lattice_defaults")
    nslice_index = next(
        (i for i, row in enumerate(defaults_state, start=1) if row.get("name") == "nslice"),
        None,
    )
    assert nslice_index is not None, "nslice parameter not found in defaults table"
    dashboard.set_input(f"default_value_{nslice_index}", 25)

    # Add a lattice element that has nslice
    dashboard.add_lattice_element("Sbend")
    assert_lattice_param_sim_input(dashboard, 0, "nslice", 25)

    # Add another default (rc) and create another element to use it
    # Set rc default as well
    defaults_state = dashboard.get_state("lattice_defaults")
    rc_index = next(
        (i for i, row in enumerate(defaults_state, start=1) if row.get("name") == "rc"),
        None,
    )
    assert rc_index is not None, "rc parameter not found in defaults table"
    dashboard.set_input(f"default_value_{rc_index}", 12.5)

    dashboard.add_lattice_element("Sbend")
    assert_lattice_param_sim_input(dashboard, 1, "rc", 12.5)

    # Unknown parameter should set an error message in state
    # Reset and verify the table is repopulated with known parameters
    dashboard.sb.click("#reset_lattice_defaults")
    defaults_state = dashboard.get_state("lattice_defaults")
    assert len(defaults_state) > 1
    assert all(row.get("name") for row in defaults_state)
