"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import importlib
import pytest
from seleniumbase import SB

from .examples import DashboardExamples

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

            # test lattice builds w/importing
            test_examples.chicane_lattice() # uses .reverse(), limited to 11 elements - ids not showing past 11
            test_examples.dogleg_lattice() # uses lattice.append() and lattice.extend()
            test_examples.expanding_fft_lattice() # uses blocking factor, and .extend() is a list of variables and direct elements
            # test_examples.iotalattice() # stress test - dashboard does not yet success well
    finally:
        if app_process is not None:
            app_process.terminate()