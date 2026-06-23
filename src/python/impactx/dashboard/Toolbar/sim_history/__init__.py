from ... import state
from .components import SimulationHistoryComponents
from .dialogs import SimulationHistoryDialogs


def save_view_details_log():
    state.sims[state.sim_index]["log"] = state.curr_view_details_log


__all__ = [
    "SimulationHistoryDialogs",
    "SimulationHistoryComponents",
]
