"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

def test_lattice_stats(dashboard):
    """
    Tests the lattice statistics to ensure accuracy and indirectly
    tests the dashboard's python parser to correctly import all inputs
    from a loaded Python file.
    """

    dashboard.load_example("fodo/run_fodo.py")

    FODO_EXPECTED_RESULTS = {
        "total_elements": 11,
        "total_length": "3.0m",
        "max_length": "1.0m",
        # "avg_length" : "1.67m",
        "min_length": "0.25m",
        "total_steps": 125,
        "periods": 1,
        "element_counts": {'beammonitor': 6, 'drift': 3, 'quad': 2}
    }


    for state_name, expected_value in FODO_EXPECTED_RESULTS.items():
        dashboard.assert_state(state_name, expected_value)

    dashboard.load_example("run/apro.py")
    APOCHROMATIC_EXPECTED_RESULTS = {
        "total_elements": 14,
        "total_length": "32.5m",
        "total_steps": 300,
        "periods": 1,
        "element_counts": {'chrquad': 8, 'chrdrift': 4, 'beammonitor': 2}
    }

    for state_name, expected_value in APOCHROMATIC_EXPECTED_RESULTS.items():
        dashboard.assert_state(state_name, expected_value)

    dashboard.load_example("chicane/run_chicane_csr.py")
    APOCHROMATIC_EXPECTED_RESULTS = {
        "total_elements": 14,
        "total_length": "32.5m",
        "total_steps": 300,
        "periods": 1,
        "element_counts": {'chrquad': 8, 'chrdrift': 4, 'beammonitor': 2}
    }

    for state_name, expected_value in APOCHROMATIC_EXPECTED_RESULTS.items():
        dashboard.assert_state(state_name, expected_value)


    dashboard.load_example("dogleg/run_dogleg_reverse.py")
    FODO_EXPECTED_RESULTS = {
        "total_elements": 11,
        "total_length": "3.0m",
        "total_steps": 125,
        "periods": 1,
        "element_counts": {'beammonitor': 6, 'drift': 43, 'quad': 2}
    }
