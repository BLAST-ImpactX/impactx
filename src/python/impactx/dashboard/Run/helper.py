import asyncio
import sys

from .. import setup_server

server, state, ctrl = setup_server()


class SimulationHelper:
    """
    Methods to help factilitate proper ImpactX simulation
    excution on the dashboard.
    """

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

        state.simulation_running = False
        ctrl.terminal_print("Simulation complete.")
        state.dirty("filtered_data")
        state.flush()


class SimulationProgress:
    """
    Methods which facilitate providing the dashboard user
    simulation progress
    """

    @staticmethod
    def print_to_xterm(content: str) -> None:
        """
        Prints the simulation content to the dashboard terminal.
        """

        ctrl.terminal_print(content.strip())
