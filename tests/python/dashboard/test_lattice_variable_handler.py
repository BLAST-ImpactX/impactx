"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import pytest


def lattice_value(state, index: int, param_name: str) -> float:
    """
    Returns the simulation value for a specific parameter
    within a lattice element at the given index.
    """

    parameters = state[index]["parameters"]
    for param in parameters:
        if param["parameter_name"] == param_name:
            value = param["sim_input"]
            if not isinstance(value, (int, float)):
                raise AssertionError(f"'{param_name}' in index {index} is not numeric.")
            return value

    raise AssertionError(f"Element {index} has no '{param_name}'")


def variable_index(variable_name: str, dashboard) -> int:
    """
    Returns the index of the variable name.
    """
    for index, variable in enumerate(dashboard.get_state("variables"), start=1):
        if variable["name"] == variable_name:
            return index
    raise ValueError(f"Variable '{variable_name}' not found in dashboard state")


def test_lattice_variable_handler(dashboard):
    """
    Verifies the lattice variable handler

    1. Add lattice element and map its variables.
    2. Delete a variable and check that the corresponding lattice parameter is invalid.
    3. Update a variable's value and ensure the lattice element reflects the change.
    4. Confirm invalid variable names produce an error.
    5. Reset variables and validate.
    """

    LATTICE_INPUTS = {"ds1": "lb", "rc1": "-rc", "nslice1": "ns"}

    VALID_VARIABLES = [
        {"name": "ns", "value": 25},
        {"name": "rc", "value": 10.3462283686195526},
        {"name": "lb", "value": 0.500194828041958},
    ]

    INVALID_VARIABLES = [
        {"name": "", "value": "3*2"},  # empty name, str
        {"name": "123abc"},  # starts with numbers
        {"name": "!"},  # special char name
        {"name": "with space"},  # space in name
        {"name": "x$"},  # symbol in name
        {"name": "rc"},  # duplicate name
    ]
    ALL_VARIABLES = VALID_VARIABLES + INVALID_VARIABLES

    # Enter lattice element and its inputs
    dashboard.add_lattice_element("Sbend")
    for element_id, value in LATTICE_INPUTS.items():
        dashboard.set_input(element_id, value)

    # Add variable rows
    dashboard.sb.click("#lattice_settings")
    for i in range(len(ALL_VARIABLES) - 1):
        dashboard.sb.click(f"#add_variable_button_{i + 1}")

    # Enter variables on the dashboard
    for i, var in enumerate(ALL_VARIABLES, start=1):
        dashboard.set_input(f"variable_name_{i}", var.get("name", ""))
        dashboard.set_input(f"variable_value_{i}", var.get("value", ""))

    # Delete 'ns' variable
    ns_var_index = variable_index("ns", dashboard)
    dashboard.sb.click(f"#delete_variable_button_{ns_var_index}")

    lattice_list = dashboard.get_state("selected_lattice_list")
    variables = dashboard.get_state("variables")

    # Verify the lattice element value is storing the variable value in the backend
    assert lattice_value(lattice_list, 0, "ds") == 0.500194828041958
    assert lattice_value(lattice_list, 0, "rc") == -10.346228368619553
    assert all(var["name"] != "ns" for var in variables), (
        "Variable 'ns' was not deleted"
    )

    # Make sure nslice parameter does not have numeric value
    # since 'ns' is deleted
    with pytest.raises(AssertionError, match="'nslice'"):
        lattice_value(lattice_list, 0, "nslice")

    # Verify that changing a variable value correctly updates the lattice element
    NEW_LB = 0.75
    lb_var_index = variable_index("lb", dashboard)
    dashboard.set_input(f"variable_value_{lb_var_index}", NEW_LB)
    lattice_list_with_updated_lb = dashboard.get_state("selected_lattice_list")
    assert lattice_value(lattice_list_with_updated_lb, 0, "ds") == NEW_LB

    # Verify that the invalid variables cannot be set on the dashboard
    for var in variables[len(VALID_VARIABLES) :]:
        assert var.get("error_message"), f"Expected an error message from {var}."

    # Reset variables and check state
    dashboard.sb.click("#reset_variables")
    variables_after_reset = dashboard.get_state("variables")
    assert len(variables_after_reset) == 1, (
        "Resetting variables should leave only one variable."
    )
    assert variables_after_reset[0].get("name", "") == "", (
        "Variable name should be blank after reset."
    )
    assert variables_after_reset[0].get("value", "") == "", (
        "Variable value should be blank after reset."
    )
    assert dashboard.get_state("is_only_variable") is True, (
        "Expected state to be True after reset."
    )
