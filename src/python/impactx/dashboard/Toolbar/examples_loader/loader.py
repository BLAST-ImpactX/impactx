"""
This file is part of ImpactX
Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from pathlib import Path

from ... import state
from ...Input.utils import GeneralFunctions
from ..file_imports.ui_populator import populate_impactx_simulation_file_to_ui

state.impactx_example_list = []


DASHBOARD_EXAMPLES = {
    "fodo/run_fodo.py",
    # "chicane/run_chicane_csr.py",
    # "fodo_space_charge/run_fodo_envelope_sc.py",
    # "apochromatic/run_apochromatic.py",
    # "kurth/run_kurth_10nC_periodic.py",
    # "expanding_beam/run_expanding_fft.py",
    # "expanding_beam/run_expanding_envelope.py",
    # "iota_lattice/run_iotalattice.py",
    # "cyclotron/run_cyclotron.py",
    # "dogleg/run_dogleg.py",
}


class DashboardExamplesLoader:
    @state.change("impactx_example")
    def on_selected_impactx_example_change(**kwargs):
        if state.impactx_example is None:
            GeneralFunctions.reset_inputs("all")
        if state.impactx_example in state.impactx_example_list:
            example_script = DashboardExamplesLoader._get_example_content(
                state.impactx_example
            )
            if example_script:
                populate_impactx_simulation_file_to_ui(example_script)

    @staticmethod
    def get_impactx_path() -> Path:
        """
        Helper method to find the impactx/examples parent directory.
        For now, just utilized to load the file names of the impactx examples
        and retrieve the script in string format.
        Potentially, this could be used to locate other impactx directories.
        """

        current_directory = Path(__file__).resolve()

        for parent in current_directory.parents:
            desired_path = parent / "examples"
            if parent.name == "impactx" and desired_path.is_dir():
                return parent

        return None

    @staticmethod
    def _get_example_content(file_name: str) -> dict:
        """
        Retrieve the selected ImpactX example file and populate the UI with its values.
        """

        impactx_directory = DashboardExamplesLoader.get_impactx_path()
        impactx_example_file_path = impactx_directory / "examples" / file_name

        file_content_as_str = impactx_example_file_path.read_text()
        file_dict = {"content": file_content_as_str.encode("utf-8")}

        return file_dict

    @staticmethod
    def load_impactx_examples() -> None:
        """
        Loads only the ImpactX example files defined in DASHBOARD_EXAMPLES
        to state.impact_example_list.
        """

        state.impactx_example_list.clear()

        impactx_directory = DashboardExamplesLoader.get_impactx_path()
        impactx_examples_directory = impactx_directory / "examples"

        for path in impactx_examples_directory.glob("**/run*"):
            relative_path = path.relative_to(impactx_examples_directory)
            relative_str = str(relative_path)
            if relative_str in DASHBOARD_EXAMPLES:
                state.impactx_example_list.append(relative_str)
