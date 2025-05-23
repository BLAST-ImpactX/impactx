"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import os
import subprocess
import sys
import time
from pathlib import Path

from selenium.common.exceptions import TimeoutException

TIMEOUT = 60


def start_dashboard() -> subprocess.Popen[str]:
    """
    Starts the impactx-dashboard in a subprocess.
    """
    repo_root = get_impactx_root_dir()
    working_directory = os.path.normpath(
        os.path.join(repo_root, "src", "python", "impactx")
    )

    print(f"Starting dashboard server from: {working_directory}")
    return subprocess.Popen(
        [sys.executable, "-u", "-m", "dashboard", "--server", "--port", "8080"],
        cwd=working_directory,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )


def get_impactx_root_dir():
    """
    Locates the ImpactX source directory.

    Looks for the first parent directory named 'impactx' that contains a '.git' folder.
    """

    current_directory = Path(__file__).resolve()

    for parent_dir in current_directory.parents:
        if parent_dir.name == "impactx" and (parent_dir / ".git").is_dir():
            return parent_dir
    return None


def wait_for_interaction_ready(sb, timeout=TIMEOUT):
    """
    Waits until the dashboard has fully loaded and is ready for interaction.

    Helper function from:
    https://github.com/Kitware/trame-client/blob/master/trame_client/utils/testing.py#L132

    """
    for i in range(timeout):
        print(f"Waiting for dashboard to load - ({i + 1}s elapsed)")
        if sb.is_element_present(".trame__loader"):
            sb.sleep(1)
        else:
            print("Ready to interact with.")
            return


def wait_for_server_ready(process, timeout=TIMEOUT):
    """
    Waits for the server of the dashboard to launch.
    """
    for i in range(timeout):
        print(f"Waiting for dashboard server to start - {i}s elapsed.")
        line = process.stdout.readline()
        if line:
            print(line, end="")
            if "App running at:" in line:
                for _ in range(2):
                    next_line = process.stdout.readline()
                    if next_line:
                        print(next_line, end="")
                print("Dashboard server has successfully launched.")
                return
        else:
            time.sleep(1)

    raise RuntimeError(f"Dashboard server did not start within {timeout} seconds.")


class DashboardTester:
    def __init__(self, sb):
        self.sb = sb

    def add_lattice_element(self, element_name: str) -> None:
        """
        Add a lattice element to the dashboard.

        :param element_name: Name of the lattice element to add.
        """
        try:
            self.set_state("selected_lattice", element_name)
            self.sb.click("#add_lattice_element")
        except Exception as error:
            raise Exception(
                f"Unable to set input for lattice element '{element_name}': {str(error)}"
            )

    def set_js_input(self, element_id: str, new_input) -> None:
        """
        Set input using JavaScript by accessing the DOM.

        Primarily used nested states, such as distribution parameters or the lattice configuration.
        However, this function can still be utilized for inputs using the 'set_state()' method.

        :param element_id: ID of the input element to set. The id is the same as the v_model_name.
        :param new_input: New value to set for the input element.
        """
        try:
            self.sb.execute_script(
                f'document.getElementById("{element_id}").value = "{new_input}";'
            )
            self.sb.execute_script(
                f'document.getElementById("{element_id}").dispatchEvent(new Event("input"));'
            )
        except Exception as error:
            raise Exception(
                f"Unable to set input for lattice element '{element_id}': {str(error)}"
            )

    def set_state(self, state_name: str, state_value):
        """
        Sets the given state name to the specified value.

        :param state_name: Name of the state variable to update (same as v_model_name).
        :param state_value: Value to set for the state.
        """
        js_script = """
            const state = window.trame?.state;
            const state_name = arguments[0];
            const state_value = arguments[1];
            if (state?.set) { state.set(state_name, state_value); }
            else if (state) {
                state[state_name] = state_value;
                if (state.dirty) state.dirty(state_name);
            }
        """
        self.sb.execute_script(js_script, state_name, state_value)

    def set_input(self, element_id: str, new_value) -> None:
        """
        Set input value by first trying set_state, then falling back to set_js_input if needed.

        This method combines the functionality of set_state and set_js_input, providing
        a unified interface for setting input values in the dashboard.

        :param element_id: ID of the input element to set (same as v_model_name).
        :param new_value: New value to set for the input element.
        """
        # Check if the state exists in trame state
        js_check_state = """
            if (window.trame && window.trame.state) {
                const state = window.trame.state;
                const state_name = arguments[0];
                if (state.get) {
                    return state.get(state_name) !== undefined;
                }
                return state.hasOwnProperty(state_name);
            }
            return false;
        """

        try:
            state_exists = self.sb.execute_script(js_check_state, element_id)
        except:
            state_exists = False

        if state_exists:
            self.set_state(element_id, new_value)
        else:
            self.set_js_input(element_id, new_value)

    def assert_state(self, state_name: str, expected_input, timeout=TIMEOUT):
        """
        Checks if trame.state[state_name] == expected_input.

        :param state_name: Name of the state variable to update (same as v_model_name).
        :param expected_input: Expected value of the state variable.
        """
        js_script = """
            if (window.trame && window.trame.state) {
                const state = window.trame.state;
                const state_name = arguments[0];
                if (state.get) { return state.get(state_name); }
                return state[state_name];
            }
            return null;
        """
        for i in range(timeout):
            try:
                value = self.sb.execute_script(js_script, state_name)
            except TimeoutException:
                value = None

            if isinstance(expected_input, (int, float)):
                if value is not None and float(value) == float(expected_input):
                    return
            elif value == expected_input:
                return

            print(
                f"Waiting for state['{state_name}'] to become '{expected_input}' - ({i + 1}s elapsed)"
            )
            time.sleep(1)

        raise TimeoutError(
            f"state['{state_name}'] never became '{expected_input}' after {timeout} seconds (last value: '{value}')."
        )
