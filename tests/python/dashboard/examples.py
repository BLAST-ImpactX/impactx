import os
from pathlib import Path
from .utils import DashboardTester, get_impactx_root_dir

class DashboardExamples:
    def __init__(self, sb):
        self.sb = sb
        self.dashboard = DashboardTester(sb)
        
        # Set up the examples directory path once
        impactx_directory = Path(get_impactx_root_dir())
        self.examples_directory = impactx_directory / "examples"
    
    def _load_example(self, example_path):
        """
        Load an example file into the dashboard
        """
        full_path = os.path.join(str(self.examples_directory), example_path)
        self.sb.click("#reset_all_inputs_button")
        self.sb.choose_file('input[type="file"]', full_path)

    def _clear(self):
        self.sb.click("#reset_all_inputs_button")

    def _begin(self, example_name: str):
        print(f"Starting {example_name} test.")

    def _end(self, example_name: str):
        print(f"Successfully ending {example_name} test.")

    def _test_lattice(self, inputs: list):
        for element_id, expected_value in inputs:
            actual_value = self.sb.get_value(element_id)

            # Try numeric comparison if expected_value is a number
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
            ("#ds2", "lb"),
            ("#rc2", "-rc"),
            ("#nslice2", "ns"),

            ("#name3", "dipedge1"),
            ("#psi3", "-psi"),
            ("#rc3", "-rc"),

            ("#name4", "dr1"),
            ("#ds4", 5.0058489435),
            ("#nslice4", "ns"),

            ("#name5", "dipedge2"),
            ("#psi5", "psi"),
            ("#rc5", "rc"),

            ("#name6", "sbend2"),
            ("#ds6", "lb"),
            ("#rc6", "rc"),
            ("#nslice6", "ns"),

            ("#name7", "dr2"),
            ("#ds7", 1.0),
            ("#nslice7", "ns"),

            # SECOND LATTICE_HALF (mirror of above, reversed, without reusing keys)
            ("#name8", "sbend2"),
            ("#ds8", "lb"),
            ("#rc8", "rc"),
            ("#nslice8", "ns"),

            ("#name9", "dipedge2"),
            ("#psi9", "psi"),
            ("#rc9", "rc"),

            ("#name10", "dr1"),
            ("#ds10", 5.0058489435),
            ("#nslice10", "ns"),

            ("#name11", "dipedge1"),
            ("#psi11", "-psi"),
            ("#rc11", "-rc"),

            ("#name12", "sbend1"),
            ("#ds12", "lb"),
            ("#rc12", "-rc"),
            ("#nslice12", "ns"),
        ]

        self._run_example("chicane", "chicane/run_chicane.py", LATTICE_CONFIGURATION)
    