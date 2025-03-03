import asyncio
import re

from .. import setup_server
from . import SimulationHelper, SimulationProgress
from .simulation import input_file

server, state, ctrl = setup_server()

start_timer = 0
state.current_step = 0

def run_execute_impactx_sim():
    asyncio.get_running_loop().create_task(execute_impactx_sim())


async def execute_impactx_sim() -> None:
    """
    Executes an ImpactX simulation based on the dashboard inputs.

    Upon call, gathers dashboard inputs, launches the simulation as
    an async subprocess, and streams its output to the dashboard terminal
    in real time.
    """

    simulation_contents = input_file()
    simulation_process = await SimulationHelper.run_simulation_in_subprocess(
        simulation_contents
    )

    while True:
        sim_output_line = await simulation_process.stdout.readline()
        sim_output_line_decoded = sim_output_line.decode()

        if not sim_output_line:
            break
        
        if "Initializing AMReX" in sim_output_line_decoded:
            start_timer = asyncio.create_task(SimulationProgress.dashboard_timer())
        if "++++ Starting step=" in sim_output_line_decoded:
            match = re.search(r"\+\+\+\+ Starting step=(\d+)", sim_output_line_decoded)
            if match:
                state.current_step = int(match.group(1))

        SimulationProgress.print_to_xterm(sim_output_line)

    await simulation_process.wait()
    start_timer.cancel()
    SimulationHelper.complete_simulation()
