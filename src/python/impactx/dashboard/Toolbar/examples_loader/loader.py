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
state.is_loading_impactx_example = False
state.has_loaded_impactx_example = False

IMPACTX_PYTHON_ROOT = Path(__file__).resolve().parents[3]

DASHBOARD_EXAMPLES = {
    "fodo/run_fodo.py",
    "chicane/run_chicane_csr.py",
    "fodo_space_charge/run_fodo_envelope_sc.py",
    "apochromatic/run_apochromatic.py",
    # "kurth/run_kurth_10nC_periodic.py", - running into recursion issues
    "expanding_beam/run_expanding_fft.py",
    "expanding_beam/run_expanding_envelope.py",
    "iota_lattice/run_iotalattice.py",
    "cyclotron/run_cyclotron.py",
    "dogleg/run_dogleg.py",
}


class DashboardExamplesLoader:
    @state.change("impactx_example")
    def on_selected_impactx_example_change(**kwargs):
        if state.is_loading_impactx_example:
            return

        selected_example = state.impactx_example

        if selected_example is None:
            if state.has_loaded_impactx_example:
                state.has_loaded_impactx_example = False
                GeneralFunctions.reset_inputs("all")
            return

        if selected_example in state.impactx_example_list:
            example_script = DashboardExamplesLoader._get_example_content(
                selected_example
            )
            if example_script:
                state.is_loading_impactx_example = True
                try:
                    populate_impactx_simulation_file_to_ui(example_script)
                    state.has_loaded_impactx_example = True
                finally:
                    state.is_loading_impactx_example = False

    @staticmethod
    def get_examples_directory() -> Path:
        """
        Return the packaged examples directory when available, otherwise fall
        back to the repo-root examples tree for source-tree development.
        """

        # in wheel-installed packages, go to examples/
        packaged_examples = IMPACTX_PYTHON_ROOT / "examples"
        if packaged_examples.is_dir():
            return packaged_examples

        # in the repo, go to two up: src/python/../../examples
        repo_examples = IMPACTX_PYTHON_ROOT.parents[2] / "examples"
        if repo_examples.is_dir():
            return repo_examples

        raise FileNotFoundError("Unable to locate ImpactX examples directory")

    @staticmethod
    def _get_example_content(file_name: str) -> dict:
        """
        Retrieve the selected ImpactX example file and populate the UI with its values.
        """

        impactx_example_file_path = (
            DashboardExamplesLoader.get_examples_directory() / file_name
        )

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

        impactx_examples_directory = DashboardExamplesLoader.get_examples_directory()

        for path in impactx_examples_directory.glob("**/run*"):
            relative_path = path.relative_to(impactx_examples_directory)
            relative_str = str(relative_path)
            if relative_str in DASHBOARD_EXAMPLES:
                state.impactx_example_list.append(relative_str)
