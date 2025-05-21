"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from datetime import datetime
from pathlib import Path

from ... import html, setup_server, vuetify
from ...Run.simulation import dashboard_sim_inputs
from . import SimulationHistoryComponents, SimulationHistoryDialogs

server, state, ctrl = setup_server()
from ..importParser import DashboardParser

state.curr_view_details_log = ""

state.sims = []
state.filtered_sims = []
state.selected_sim_to_load = None
state.sim_to_download = None
state.selected_sim_search_status = "All"

state.sim_history_table_headers = [
    {"title": "Simulation Name", "key": "name", "sortable": True},
    {"title": "Created", "key": "created_at_time", "sortable": True, "align": "center"},
    {"title": "Duration", "key": "time_elapsed", "sortable": True, "align": "center"},
    {"title": "Status", "key": "status", "sortable": True, "align": "center"},
    {"title": "Actions", "key": "actions", "sortable": False, "align": "center"},
]


# --------------------------------
# Load custom JS
# --------------------------------


def load_my_js(server):
    """
    Loads custom js file to the server.
    """
    js_file = Path(__file__).with_name("custom.js").resolve()
    server.enable_module(
        {
            "serve": {"my_code": str(js_file.parent)},
            "scripts": [f"my_code/{js_file.name}"],
        }
    )


# --------------------------------
# Functionality
# --------------------------------


class SimulationHistory:
    """
    Builds the UI and handles functionality to handle
    simulation history for the dashboard.
    """

    @staticmethod
    @ctrl.add("open_view_details")
    def open_view_details(selected_sim):
        state.selected_sim = selected_sim
        state.view_details_dialog = True

    @staticmethod
    @ctrl.add("rename_sim")
    def open_rename_dialog(sim):
        state.rename_old_name = sim["name"]
        state.rename_new_name = sim["name"]
        state.sim_rename_dialog = True

    @staticmethod
    @ctrl.add("confirm_rename")
    def confirm_rename():
        old = state.rename_old_name
        new = state.rename_new_name

        for sim in state.sims:
            if sim["name"] == old:
                sim["name"] = new
                break

        state.rename_dialog = False
        state.dirty("filtered_sims")

        SimulationHistory._close_rename_dialog()

    @staticmethod
    @ctrl.add("close_rename_dialog")
    def _close_rename_dialog():
        state.sim_rename_dialog = False

    @staticmethod
    @ctrl.add("update_search")
    def update_search(user_input):
        state.selected_sim_search = user_input
        SimulationHistory.filter_sim_history()

    @staticmethod
    @ctrl.add("update_status")
    def update_status(user_input):
        state.selected_sim_status = user_input
        SimulationHistory.filter_sim_history()

    @staticmethod
    @ctrl.trigger("download_sim")
    def download_sim(sim):
        return sim.get("inputs", "")

    @staticmethod
    @ctrl.add("delete_sim")
    def delete_sim(sim):
        actual_index = state.sims.index(sim)

        del state.sims[actual_index]
        SimulationHistory.filter_sim_history()
        state.dirty("filtered_sims")
        state.dirty("sims")

    @staticmethod
    @ctrl.add("toggle_selected_sim")
    def toggle_selected_sim(sim):
        same_name = (
            state.selected_sim_to_load is not None
            and state.selected_sim_to_load["name"] == sim["name"]
        )

        if same_name:
            state.selected_sim_to_load = None
        else:
            state.selected_sim_to_load = sim

    @staticmethod
    @ctrl.add("load_selected_sim")
    def load_selected_sim():
        sim = state.selected_sim_to_load
        sim_content = {
            "name": sim["name"] + ".py",
            "content": sim["inputs"].encode("utf-8"),
        }

        DashboardParser.populate_impactx_simulation_file_to_ui(sim_content)
        state.selected_sim_to_load = None

    @staticmethod
    def filter_sim_history():
        """
        Handles the functionality to filter the sim history.
        """
        filtered = state.sims

        if (
            state.selected_sim_search_status
            and state.selected_sim_search_status != "All"
        ):
            filtered = [
                sim
                for sim in filtered
                if sim["status"] == state.selected_sim_search_status
            ]

        if state.selected_sim_search:
            search_query = state.selected_sim_search.lower()
            filtered = [sim for sim in filtered if search_query in sim["name"].lower()]

        state.filtered_sims = filtered

    # --------------------------------
    # Helper Functions
    # --------------------------------

    @staticmethod
    def _access_sim_history_slot(key):
        return vuetify.Template(raw_attrs=[f'v-slot:item.{key}="{{ item }}"'])

    @staticmethod
    def init_sim_history_dialogs():
        SimulationHistoryDialogs.rename_dialog()
        SimulationHistoryDialogs.view_details_dialog()
        SimulationHistoryDialogs.download_options_dialog()
        SimulationHistoryDialogs.load_sim_dialog()

    @staticmethod
    def add_sim_to_history():
        """
        Called once a simulation is ran.

        Adds a new simulation to the sim history.
        """

        curr_num_sims = len(state.sims)
        new_sim_name = state.imported_file_name or f"Simulation_{curr_num_sims + 1}"
        current_time = datetime.utcnow().isoformat() + "Z"
        sim_inputs = dashboard_sim_inputs(is_exporting=True)

        new_sim = {
            "name": new_sim_name,
            "created_at_time": current_time,
            "time_elapsed": "",
            "status": "In Progress",
            "inputs": sim_inputs,
            "log": "",
        }

        state.sims = state.sims + [new_sim]
        state.filtered_sims = state.sims
        return curr_num_sims

    @staticmethod
    def add_to_view_details_log(log: str) -> None:
        """
        Stores simulation details inside of
        curr_view_details_log state
        """

        state.curr_view_details_log += log

    @staticmethod
    def save_view_details_log() -> None:
        """
        Updates the UI with the simulation's
        details.
        """

        state.sims[state.sim_index]["log"] = state.curr_view_details_log

    @staticmethod
    def simulation_history():
        """
        Contains the UI and functionality for the
        dashboard's simulation history.
        """
        with SimulationHistoryComponents.dialog(
            title="Simulation History",
            prepend_icon="mdi-clipboard-text-clock",
            dialog_var="simulation_history_dialog",
            width="75vw",
        ):
            with vuetify.VCardText():
                with vuetify.VRow():
                    with vuetify.VCol(cols=12, sm=8):
                        SimulationHistoryComponents.text_field(
                            label="Search simulations",
                            v_model_name="selected_sim_search",
                            update_modelValue=(ctrl.update_status, "[$event]"),
                            prepend_inner_icon="mdi-magnify",
                            clearable=True,
                        )
                    with vuetify.VCol(cols=12, sm=4):
                        vuetify.VSelect(
                            label="Status",
                            v_model=("selected_sim_search_status", None),
                            update_modelValue=(ctrl.update_status, "[$event]"),
                            items=(
                                [
                                    "All",
                                    "Completed",
                                    "In Progress",
                                    "Cancelled",
                                    "Failed",
                                ],
                            ),
                            clearable=True,
                            density="comfortable",
                            hide_details=True,
                            variant="outlined",
                        )
                with vuetify.VRow():
                    with vuetify.VCol(cols=12):
                        with vuetify.VDataTable(
                            classes="elevation-2",
                            headers=("sim_history_table_headers",),
                            items=("filtered_sims",),
                        ):
                            with SimulationHistory._access_sim_history_slot("name"):
                                with html.Div(
                                    style="display: flex; align-items: center; gap: 6px;"
                                ):
                                    html.Span("{{ item.name }}")
                                    SimulationHistoryComponents.icon_button(
                                        icon_name="mdi-pencil",
                                        color="warning",
                                        click=(ctrl.rename_sim, "[item]"),
                                        description="Rename",
                                    )
                            with SimulationHistory._access_sim_history_slot(
                                "created_at_time"
                            ):
                                with html.Div(
                                    style="display: flex; flex-direction: column; align-items: center;"
                                ):
                                    html.Div(
                                        "{{ window.formatDate(item.created_at_time) }}",
                                        style="font-weight: 500",
                                    )
                                    html.Div(
                                        "{{ window.formatTime(item.created_at_time) }}",
                                        classes="text-caption",
                                    )
                            with SimulationHistory._access_sim_history_slot(
                                "time_elapsed"
                            ):
                                html.Span("{{ item.time_elapsed || '—' }}")
                            with SimulationHistory._access_sim_history_slot("status"):
                                SimulationHistoryComponents.status_chip("item")
                            with SimulationHistory._access_sim_history_slot("actions"):
                                SimulationHistoryComponents.icon_button(
                                    icon_name="mdi-eye",
                                    classes="mr-1",
                                    click=(ctrl.open_view_details, "[item]"),
                                    description="View Details",
                                    disabled=("sim_is_running",),
                                )
                                SimulationHistoryComponents.icon_button(
                                    icon_name="mdi-download",
                                    click="""
                                            sim_to_download = item;
                                            sim_download_dialog = true;
                                        """,
                                    description="Download",
                                )
                                SimulationHistoryComponents.icon_button(
                                    icon_name="mdi-tray-arrow-up",
                                    color="primary",
                                    description="Load",
                                    click="""
                                                selected_sim_to_load = item;
                                                load_sim_dialog = true;
                                            """,
                                )
                                SimulationHistoryComponents.icon_button(
                                    icon_name="mdi-trash-can-outline",
                                    color="error",
                                    click=(ctrl.delete_sim, "[item]"),
                                    description="Delete",
                                    disabled=("sim_is_running",),
                                )
