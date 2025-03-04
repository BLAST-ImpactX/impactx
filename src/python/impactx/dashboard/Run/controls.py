import asyncio
import re

from .. import setup_server
from . import SimulationHelper, SimulationProgress
from .simulation import input_file

server, state, ctrl = setup_server()

state.sim_elapsed_time = "0.0"
state.sim_is_running = False
state.sim_current_step = 0
state.sim_total_steps = 0
state.sim_progress = 0
start_timer = 0

def run_execute_impactx_sim():
    asyncio.get_running_loop().create_task(execute_impactx_sim())
state.sim_status_color = "primary"


async def execute_impactx_sim() -> None:
    """
    Executes an ImpactX simulation based on the dashboard inputs.

    Upon call, gathers dashboard inputs, launches the simulation as
    an async subprocess, and streams its output to the dashboard terminal
    in real time.
    """

    simulation_contents = input_file()
    state.sim_total_steps = SimulationProgress.determine_sim_total_steps(simulation_contents)
    simulation_process = await SimulationHelper.run_simulation_in_subprocess(
        simulation_contents
    )

    while True:
        state.sim_is_running = True
        sim_output_line = await simulation_process.stdout.readline()
        sim_output_line_decoded = sim_output_line.decode()

        if not sim_output_line:
            break
        
        if "Initializing AMReX" in sim_output_line_decoded:
            start_timer = asyncio.create_task(SimulationProgress.dashboard_timer())
        if "++++ Starting step=" in sim_output_line_decoded:
            match = re.search(r"\+\+\+\+ Starting step=(\d+)", sim_output_line_decoded)
            if match:
                state.sim_current_step = int(match.group(1))
                state.sim_progress = (state.sim_current_step / state.sim_total_steps) * 100

        SimulationProgress.print_to_xterm(sim_output_line)

    await simulation_process.wait()
    start_timer.cancel()
    SimulationHelper.complete_simulation()
