from ... import setup_server
from .components import SimulationHistoryComponents
from .dialogs import SimulationHistoryDialogs

server, state, ctrl = setup_server()


def save_view_details_log():
    state.sims[state.sim_index]["log"] = state.curr_view_details_log


__all__ = [
    "SimulationHistoryDialogs",
    "SimulationHistoryComponents",
]
