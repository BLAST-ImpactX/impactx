import asyncio
import base64
import os
import re
import sys
import time

from .. import setup_server
from ..Toolbar.sim_history.ui import SimulationHistory

server, state, ctrl = setup_server()

state.sim_progress_status = ""


class SimulationHelper:
    """
    Methods to help factilitate proper ImpactX simulation
    excution on the dashboard.
    """

    @staticmethod
    async def run_simulation_in_subprocess(simulation_contents):
        """
        Runs the simulation script as an asynchronous subprocess.
        """

        process = await asyncio.create_subprocess_exec(
            sys.executable,
            "-c",
            simulation_contents,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT,
        )

        return process

    @staticmethod
    def complete_simulation():
        """
        Marks the simulation as complete and updates the dashboard.
        """

        state.sim_is_running = False
        state.sim_progress = 100
        ctrl.terminal_print("Simulation complete.")
        state.dirty("filtered_data")
        state.sim_status_color = "success"
        state.sims[state.sim_index]["status"] = "Completed"
        state.sim_is_generating_plots = False
        state.dirty("filtered_sims")
        state.flush()
        SimulationHistory.save_view_details_log()

    @staticmethod
    def cancel_simulation(proc: asyncio.subprocess.Process):
        if proc is not None and proc.returncode is None:
            proc.kill()

        state.sim_is_cancelled = True
        state.sim_is_running = False
        state.sim_progress = 0
        state.sim_current_step = 0
        state.sim_elapsed_time = "0.0"
        state.sim_status_color = "warning"
        state.sim_progress_status = "Cancelled"
        state.sims[state.sim_index]["status"] = "Cancelled"
        state.dirty("filtered_sims")
        ctrl.terminal_print("Simulation cancelled.")
        state.flush()
        SimulationHistory.save_view_details_log()

    @staticmethod
    def fail_simulation() -> None:
        """
        Updates the UI and simulation history to reflect a failed simulation.
        """
        state.sim_is_running = False
        state.sim_progress_status = "Failed"
        state.sim_status_color = "error"
        state.sims[state.sim_index]["status"] = "Failed"
        state.dirty("filtered_sims")
        state.flush()
        ctrl.terminal_print("Simulation failed due to the above error.")
        SimulationHistory.save_view_details_log()

    @staticmethod
    def reset():
        state.sim_is_cancelled = False
        state.sim_is_running = True
        state.sim_progress = 0
        state.sim_current_step = 0
        state.sim_elapsed_time = "0.0"
        state.sim_status_color = "primary"
        state.curr_view_details_log = ""
        state.flush()

    @staticmethod
    def display_phase_space_plots():
        """
        Loads the ImpactX generated phase space plots
        to the dashboard.
        """

        if os.path.exists("phase_space_plot.png"):
            with open("phase_space_plot.png", "rb") as f:
                image_data = f.read()
                image_base64 = base64.b64encode(image_data).decode()
                state.phase_space_png = f"data:image/png;base64, {image_base64}"
                state.flush()

            os.remove("phase_space_plot.png")


class SimulationProgress:
    """
    Methods which facilitate providing the dashboard user
    simulation progress
    """

    @state.change(
        "sim_current_step", "sim_total_steps", "sim_is_running", "sim_progress"
    )
    def _update_status(**kwargs):
        if state.sim_is_running:
            if state.sim_current_step == 0:
                state.sim_progress_status = "Starting..."
            elif state.sim_current_step >= state.sim_total_steps:
                state.sim_is_generating_plots = True
                state.sim_progress_status = "Generating plots..."
            else:
                progress_percent = int(state.sim_progress)
                state.sim_progress_status = f"Running... ({progress_percent}%)"
        elif state.sim_progress == 100:
            state.sim_progress_status = "Complete!"

    @staticmethod
    def print_to_xterm(content: str) -> None:
        """
        Prints the simulation content to the dashboard terminal.
        """

        ctrl.terminal_print(content.strip())

    @staticmethod
    async def dashboard_timer():
        start_time = time.monotonic()

        while True:
            elapsed = time.monotonic() - start_time
            state.sim_elapsed_time = SimulationProgress.format_elapsed_time(elapsed)
            state.sims[state.sim_index]["time_elapsed"] = state.sim_elapsed_time
            state.dirty("filtered_sims")
            state.flush()
            await asyncio.sleep(0.1)

    @staticmethod
    def determine_sim_total_steps(simulation_content_file) -> int:
        """
        Determines the total step count for the given input file.

        Sum of nslices is sim_total_step
        """

        nslice_matches = re.findall(r"nslice=(\d+)", simulation_content_file)

        if nslice_matches:
            state.sim_total_steps = sum(int(match) for match in nslice_matches)

        return state.sim_total_steps

    @staticmethod
    def format_elapsed_time(seconds: float) -> str:
        """
        Converts elapsed seconds to a clearly-readable string.
        """

        seconds = round(seconds, 1)
        minutes, sec = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        if hours:
            return f"{int(hours)}h {int(minutes)}m {int(sec)}s"
        elif minutes:
            return f"{int(minutes)}m {int(sec)}s"
        else:
            return f"{sec:.1f}s"
