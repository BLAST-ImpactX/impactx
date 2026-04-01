"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import time

from .utils import TIMEOUT, DashboardTester


class DashboardExamples:
    def __init__(self, dashboard: DashboardTester):
        self.dashboard = dashboard
        self.sb = dashboard.sb

    def _load_example(self, example_path):
        self.dashboard.load_example(example_path)

    def _clear(self):
        self.sb.click("#reset_all_inputs_button")

    def _begin(self, example_name: str):
        print(f"Starting {example_name} test.")

    def _end(self, example_name: str):
        print(f"Successfully ending {example_name} test.")

    def _test_lattice(self, inputs: list):
        # Use get_attribute (does not require element visibility, unlike
        # get_value which fails for scrollable containers) with a retry
        # loop to wait for async file-import state propagation.
        for element_id, expected_value in inputs:
            for _ in range(TIMEOUT):
                actual_value = self.sb.get_attribute(element_id, "value")
                if actual_value is not None and actual_value != "":
                    break
                time.sleep(1)
            else:
                raise AssertionError(
                    f"{element_id}: value was still empty after {TIMEOUT} seconds"
                )

            if isinstance(expected_value, (int, float)):
                actual_value = float(actual_value)
                assert actual_value == expected_value, (
                    f"{element_id}: expected {expected_value}, got {actual_value}"
                )
            else:
                assert actual_value == expected_value, (
                    f"{element_id}: expected '{expected_value}', got '{actual_value}'"
                )

    def _run_example(self, name: str, path: str, test_data: list):
        self._begin(name)
        self._clear()
        self._load_example(path)
        self._test_lattice(test_data)
        self._end(name)

    def dogleg_lattice(self):
        LATTICE_CONFIGURATION = [
            # sbend1
            ("#name2", "sbend1"),
            ("#ds2", "lb"),
            ("#rc2", "-rc"),
            ("#nslice2", "ns"),
            # dipedge1
            ("#name3", "dipedge1"),
            ("#psi3", "-psi"),
            ("#rc3", "-rc"),
            # dr1
            ("#name4", "dr1"),
            ("#ds4", 5.0058489435),
            ("#nslice4", "ns"),
            # dipedge2
            ("#name5", "dipedge2"),
            ("#psi5", "psi"),
            ("#rc5", "rc"),
            # sbend2
            ("#name6", "sbend2"),
            ("#ds6", "lb"),
            ("#rc6", "rc"),
            ("#nslice6", "ns"),
            # dr2
            ("#name7", "dr2"),
            ("#ds7", 0.5),
            ("#nslice7", "ns"),
        ]
        self._run_example("dogleg", "dogleg/run_dogleg.py", LATTICE_CONFIGURATION)

    def chicane_lattice(self):
        """
        Utilizes lattice_half.reverse()
        """
        LATTICE_CONFIGURATION = [
            # FIRST LATTICE_HALF
            ("#name2", "sbend1"),
            ("#name3", "dipedge1"),
            ("#name5", "dr1"),
            ("#name7", "dipedge2"),
            ("#name8", "sbend2"),
            ("#name10", "dr2"),
            ("#name12", "sbend2"),
            ("#name13", "dipedge2"),
            ("#name15", "dr1"),
            ("#name17", "dipedge1"),
            ("#name18", "sbend1"),
            ("#name20", "dr3"),
        ]

        self._run_example("chicane", "chicane/run_chicane.py", LATTICE_CONFIGURATION)


def test_examples(dashboard) -> None:
    """
    Exercise a subset of example imports against the running dashboard.
    """

    test_examples = DashboardExamples(dashboard)

    # test lattice builds w/example importing
    test_examples.dogleg_lattice()  # uses lattice.append() and lattice.extend()
    test_examples.chicane_lattice()  # uses .reverse()
