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

    def cyclotron_lattice(self):
        LATTICE_CONFIGURATION = [
            ("#name2", "acc1"),
            ("#ds2", 0.038),
            ("#ez2", 1.12188308693e-4),
            ("#bz2", 1.0e-14),
            ("#nslice2", "ns"),
            # First ExactSbend element (sbend1)
            ("#name3", "sbend1"),
            ("#ds3", 0.25),
            ("#phi3", 180.0),
            ("#B3", 1),
            # Second ChrAcc element (acc2)
            ("#name4", "acc2"),
            ("#ds4", 0.038),
            ("#ez4", 1.12188308693e-4),
            ("#bz4", 1.0e-14),
            ("#nslice4", "ns"),
            # Second ExactSbend element (sbend2)
            ("#name5", "sbend2"),
            ("#ds5", 0.25),
            ("#phi5", 180.0),
            ("#B5", 1),
        ]
        self._run_example("cyclotron", "cyclotron/run_cyclotron.py", LATTICE_CONFIGURATION)


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

    def expanding_fft_lattice(self):
        """
        Utilizes unique build
        """
        LATTICE_CONFIGURATION = [
            ("#ds2", 6.0),
            ("#nslice2", 40),
        ]
        self._run_example(
            name="expanding_fft",
            path="expanding_beam/run_expanding_fft.py",
            test_data=LATTICE_CONFIGURATION
        )

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

            # ("#name12", "sbend1"),
            # ("#ds12", "lb"),
            # ("#rc12", "-rc"),
            # ("#nslice12", "ns"),
        ]

        self._run_example("chicane", "chicane/run_chicane.py", LATTICE_CONFIGURATION)
        
    def apochromatic_lattice(self):
        LATTICE_CONFIGURATION = [
            # Drift elements
            ("#name2", "dr1"),
            ("#ds2", 1.0),
            ("#nslice2", "ns"),
            
            ("#name6", "dr2"),
            ("#ds6", 10.0),
            ("#nslice6", "ns"),
            # Quad elements
            ("#name3", "q1"),
            ("#ds3", 1.2258333333),
            ("#k3", 0.5884),
            ("#nslice3", "ns"),
            ("#name4", "q2"),
            ("#ds4", 1.5677083333),
            ("#k4", -0.7525),
            ("#nslice4", "ns"),
            ("#name5", "q3"),
            ("#ds5", 1.205625),
            ("#k5", 0.5787),
            ("#nslice5", "ns"),
            ("#name7", "q4"),
            ("#ds7", 1.2502083333),
            ("#k7", -0.6001),
            ("#nslice7", "ns"),
            ("#name8", "q5"),
            ("#ds8", 1.2502083333),
            ("#k8", 0.6001),
            ("#nslice8", "ns"),
            ("#name10", "q6"),
            ("#ds10", 1.205625),
            ("#k10", -0.5787),
            ("#nslice10", "ns"),
            
            ("#name11", "q7"),
            ("#ds11", 1.5677083333),
            ("#k11", 0.7525),
            ("#nslice11", "ns"),
            
            # ("#name12", "q8"),
            # ("#ds12", 1.2258333333),
            # ("#k12", -0.5884),
            # ("#nslice12", "ns"),
        ]
        self._run_example("apochromatic", "apochromatic/run_apochromatic.py", LATTICE_CONFIGURATION)

    def kurth_10nC_periodic_lattice(self):
        example_name = "kurth_10nC_periodic"
        LATTICE_CONFIGURATION = [
            ("#name2", "drift1"),
            ("#ds2", 1),
            ("#name3", "constf1"),
            ("#ds3", 2),
            ("#kx3", 0.7),
            ("#ky3", 0.7),
            ("#kt3", 0.7),
            ("#name4", "drift1"),
            ("#ds4", 1),
        ]
        self._run_example(
            name=example_name,
            path=f"kurth/run_{example_name}.py",
            test_data=LATTICE_CONFIGURATION
        )
        