from .components import SimulationHistoryComponents
from .dialogs import SimulationHistoryDialogs

from ... import setup_server
server, state, ctrl = setup_server()
def save_view_details_log() :
    state.sims[state.sim_index]["log"] = state.curr_view_details_log


__all__ = [
    "SimulationHistoryDialogs",
    "SimulationHistoryComponents",
]
