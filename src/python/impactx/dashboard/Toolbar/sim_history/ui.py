from pathlib import Path
from datetime import datetime

from ... import setup_server, vuetify

server, state, ctrl = setup_server()


state.simulation_history_dialog = False
state.sims = []

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
    def add_sim_to_history():
        curr_num_sims = len(state.sims)
        new_sim_name = f"Simulation_{curr_num_sims + 1}"
        current_time =  datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")

        new_sim = {
            "name": new_sim_name,
            "created_at_time": current_time,
            "time_elapsed": "2h 30m",
            "status": "In Progress",
        }

        state.sims = state.sims + [new_sim]
        return curr_num_sims

    @staticmethod
    def simulation_history():
        with vuetify.VDialog(
            v_model=("simulation_history_dialog",),
            style="width:75vw",
        ):
            with vuetify.VCard(
                elevation=8,
                classes="rounded-lg",
            ):
                with vuetify.VToolbar(
                    classes="px-4",
                    color="primary",
                ):
                    vuetify.VIcon("mdi-clipboard-text-clock")
                    vuetify.VToolbarTitle("Simulation History")
                    vuetify.VSpacer()
                    vuetify.VBtn(
                        icon="mdi-close",
                        click="simulation_history_dialog = false",
                    )
                with vuetify.VCardText():
                    with vuetify.VRow():
                        with vuetify.VCol(cols=12, sm=8):
                            vuetify.VTextField(
                                label="Search simulations",
                                prepend_inner_icon="mdi-magnify",
                                clearable=True,
                                density="comfortable",
                                hide_details=True,
                                variant="outlined",
                            )
                        with vuetify.VCol(cols=12, sm=4):
                            vuetify.VSelect(
                                label="Status",
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
                                items=("sims",),
                            ):
                                with vuetify.Template(raw_attrs=['v-slot:item.status="{ item }"']):
                                    vuetify.VChip(
                                        "{{ item.status }}",
                                        color=("window.getSimStatusColor(item.status)",),
                                        variant="elevated",
                                    )
