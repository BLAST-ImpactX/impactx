"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import importlib
import pytest
from seleniumbase import SB

from examples import DashboardExamples
from utils import (
    start_dashboard,
    wait_for_interaction_ready,
    wait_for_server_ready,
)

@pytest.mark.skipif(
    importlib.util.find_spec("seleniumbase") is None,
    reason="seleniumbase is not available",
)

def test_examples():
    # Note: for the test_examples, if any of the files used for the examples change in the future, there is possibility that the example will fail
    # This isnt the worst issue as it will mean that the dahsboard will need to update accordingly.

    app_process = None
    try:
        with SB(headless=True) as sb:
            test_examples = DashboardExamples(sb)

            # Setup Dashboard
            app_process = start_dashboard()
            wait_for_server_ready(app_process)
            sb.open("http://localhost:8080/index.html#/Input")
            wait_for_interaction_ready(sb)

            # Run tests once dashboard is ready for interaction
            # Test full example inputs
            test_examples.fodo_example()

    finally:
        if app_process is not None:
            app_process.terminate()