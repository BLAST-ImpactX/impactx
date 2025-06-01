import os
from pathlib import Path
from utils import DashboardTester, get_impactx_root_dir

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
        self.sb.click("#reset_button")
        self.sb.choose_file('input[type="file"]', full_path)

    def _clear(self):
        self.sb.click("#reset_button")

    def _begin(self, example_name: str):
        print(f"Starting {example_name} test.")

    def _end(self, example_name: str):
        print(f"Successfully ending {example_name} test.")

    def _test_lattice(self, inputs: list):
        for element_id, expected_value in inputs:
            # Check if element is accessible
            # if not self._element_exists(element_id):
            #     print(f"Skipping {element_id} - not currently accessible")
            #     continue

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

    def fodo_example(self):
        self._clear()
        self._load_example("fodo/run_fodo.py")

        # Define test parameters
        BEAM_PARAMETERS = [
            ("tracking_mode", "Particle Tracking"),
            ("space_charge", "false"),
            ("charge_qe", -1),
            ("mass_MeV", 0.510998950),
            ("npart", 10000),
            ("bunch_charge_C", 1e-9)
        ]
        
        DISTRIBUTION_PARAMETERS = [
            ("distribution", "Waterbag"),
            ("distribution_type", "Quadratic"),
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
            ("#mutpt", 0.0)
        ]
        
        LATTICE_CONFIGURATION = [
            ("#ds2", 0.25),
            ("#ds4", 1),
            ("#k4", 1.0),
            ("#ds6", 0.5),
            ("#ds8", 1.0),
            ("#k8", -1.0),
            ("#ds10", 0.25)
        ]

        # Check state parameters
        for param_name, expected_value in BEAM_PARAMETERS + DISTRIBUTION_PARAMETERS:
            self.dashboard.assert_state(param_name, expected_value)

        # Check form input values
        for element_id, expected_value in DISTRIBUTION_VALUES + LATTICE_CONFIGURATION:
            actual_value = float(self.sb.get_value(element_id))
            assert actual_value == expected_value, f"{element_id}: expected {expected_value}, got {actual_value}"

