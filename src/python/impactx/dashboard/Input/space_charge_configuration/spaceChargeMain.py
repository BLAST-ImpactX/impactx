from trame.widgets import vuetify

from ...trame_setup import setup_server
from ..generalFunctions import generalFunctions
from .spaceChargeFunctions import SpaceChargeFunctions

server, state, ctrl = setup_server()

# -----------------------------------------------------------------------------
# Default
# -----------------------------------------------------------------------------

state.dynamic_size = generalFunctions.get_default("dynamic_size", "default_values")
state.max_level = generalFunctions.get_default("max_level", "default_values")
state.particle_shape = generalFunctions.get_default("particle_shape", "default_values")
state.poisson_solver = generalFunctions.get_default("poisson_solver", "default_values")

state.prob_relative = []
state.prob_relative_fields = []

state.n_cell = []
state.n_cell_x = generalFunctions.get_default("n_cell", "default_values")
state.n_cell_y = generalFunctions.get_default("n_cell", "default_values")
state.n_cell_z = generalFunctions.get_default("n_cell", "default_values")

state.blocking_factor_x = generalFunctions.get_default(
    "blocking_factor", "default_values"
)
state.blocking_factor_y = generalFunctions.get_default(
    "blocking_factor", "default_values"
)
state.blocking_factor_z = generalFunctions.get_default(
    "blocking_factor", "default_values"
)

state.mlmg_relative_tolerance = generalFunctions.get_default(
    "mlmg_relative_tolerance", "default_values"
)
state.mlmg_absolute_tolerance = generalFunctions.get_default(
    "mlmg_absolute_tolerance", "default_values"
)
state.mlmg_max_iters = generalFunctions.get_default("mlmg_max_iters", "default_values")
state.mlmg_verbosity = generalFunctions.get_default("mlmg_verbosity", "default_values")

state.error_message_mlmg_relative_tolerance = ""
state.error_message_mlmg_absolute_tolerance = ""
state.error_message_mlmg_max_iters = ""
state.error_message_mlmg_verbosity = ""

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------


def populate_prob_relative_fields(max_level):
    num_prob_relative_fields = int(max_level) + 1
    fft_first_field_value = generalFunctions.get_default(
        "prob_relative_first_value_fft", "default_values"
    )
    multigrid_first_field_value = generalFunctions.get_default(
        "prob_relative_first_value_multigrid", "default_values"
    )

    if state.poisson_solver == "fft":
        state.prob_relative = [fft_first_field_value] + [0.0] * (
            num_prob_relative_fields - 1
        )
    elif state.poisson_solver == "multigrid":
        state.prob_relative = [multigrid_first_field_value] + [0.0] * (
            num_prob_relative_fields - 1
        )
    else:
        state.prob_relative = [0.0] * num_prob_relative_fields

    state.prob_relative_fields = [
        {
            "value": state.prob_relative[i],
            "error_message": SpaceChargeFunctions.validate_prob_relative_fields(
                i, state.prob_relative[i]
            ),
            "step": generalFunctions.get_default("prob_relative", "steps"),
        }
        for i in range(num_prob_relative_fields)
    ]


def update_blocking_factor_and_n_cell(category, kwargs):
    directions = ["x", "y", "z"]

    for state_name, value in kwargs.items():
        if any(state_name == f"{category}_{dir}" for dir in directions):
            direction = state_name.split("_")[-1]
            SpaceChargeFunctions.validate_n_cell_and_blocking_factor(direction)

            n_cell_error = getattr(state, f"error_message_n_cell_{direction}")
            blocking_factor_error = getattr(
                state, f"error_message_blocking_factor_{direction}"
            )

            if not n_cell_error:
                n_cell_value = getattr(state, f"n_cell_{direction}")
                setattr(state, f"n_cell_{direction}", int(n_cell_value))

            if not blocking_factor_error:
                blocking_factor_value = getattr(state, f"blocking_factor_{direction}")
                setattr(
                    state, f"blocking_factor_{direction}", int(blocking_factor_value)
                )

    state.n_cell = [getattr(state, f"n_cell_{dir}", 0) for dir in directions]
    state.blocking_factor = [
        getattr(state, f"blocking_factor_{dir}", 0) for dir in directions
    ]


# -----------------------------------------------------------------------------
# Decorators
# -----------------------------------------------------------------------------
@state.change("poisson_solver")
def on_poisson_solver_change(poisson_solver, **kwargs):
    populate_prob_relative_fields(state.max_level)
    state.dirty("prob_relative_fields")
    generalFunctions.update_simulation_validation_status()


@state.change("space_charge")
def on_space_charge_change(space_charge, **kwargs):
    state.dynamic_size = space_charge
    generalFunctions.update_simulation_validation_status()


@state.change("max_level")
def on_max_level_change(max_level, **kwargs):
    populate_prob_relative_fields(max_level)
    generalFunctions.update_simulation_validation_status()


@state.change("blocking_factor_x", "blocking_factor_y", "blocking_factor_z")
def on_blocking_factor_change(**kwargs):
    update_blocking_factor_and_n_cell("blocking_factor", kwargs)
    generalFunctions.update_simulation_validation_status()


@state.change("n_cell_x", "n_cell_y", "n_cell_z")
def on_n_cell_change(**kwargs):
    update_blocking_factor_and_n_cell("n_cell", kwargs)
    generalFunctions.update_simulation_validation_status()


@ctrl.add("update_prob_relative")
def on_update_prob_relative_call(index, value):
    index = int(index)

    try:
        prob_relative_value = float(value)
        state.prob_relative[index] = prob_relative_value
    except ValueError:
        prob_relative_value = 0.0
        state.prob_relative[index] = prob_relative_value

    # Validate the updated value
    error_message = SpaceChargeFunctions.validate_prob_relative_fields(
        index, prob_relative_value
    )

    state.prob_relative_fields[index]["value"] = value
    state.prob_relative_fields[index]["error_message"] = error_message

    # Validate next index if it exists
    if index + 1 < len(state.prob_relative):
        next_value = state.prob_relative[index + 1]

        next_error_message = SpaceChargeFunctions.validate_prob_relative_fields(
            index + 1, next_value
        )
        state.prob_relative_fields[index + 1]["error_message"] = next_error_message

    state.dirty("prob_relative_fields")
    generalFunctions.update_simulation_validation_status()


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class SpaceChargeConfiguration:
    @staticmethod
    def card():
        """
        Creates UI content for space charge configuration
        """

        with vuetify.VDialog(v_model=("showSpaceChargeDialog", False), width="500px"):
            SpaceChargeConfiguration.dialog_space_charge_settings()

        with vuetify.VCard(v_show="space_charge", style="width: 340px;"):
            with vuetify.VCardTitle("Space Charge"):
                vuetify.VSpacer()
                vuetify.VIcon(
                    "mdi-cog",
                    classes="ml-2",
                    v_if="poisson_solver == 'multigrid'",
                    click="showSpaceChargeDialog = true",
                    style="cursor: pointer;",
                )
                vuetify.VIcon(
                    "mdi-information",
                    classes="ml-2",
                    click=lambda: generalFunctions.documentation(
                        "space_charge_documentation"
                    ),
                    style="color: #00313C;",
                )
            vuetify.VDivider()
            with vuetify.VCardText():
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(cols=5, classes="py-0"):
                        vuetify.VSelect(
                            label="Poisson Solver",
                            v_model=("poisson_solver",),
                            items=(
                                generalFunctions.get_default(
                                    "poisson_solver_list", "default_values"
                                ),
                            ),
                            dense=True,
                            hide_details=True,
                        )
                    with vuetify.VCol(cols=4, classes="py-0"):
                        vuetify.VSelect(
                            label="Particle Shape",
                            v_model=("particle_shape",),
                            items=(
                                generalFunctions.get_default(
                                    "particle_shape_list", "default_values"
                                ),
                            ),
                            dense=True,
                        )
                    with vuetify.VCol(cols=3, classes="py-0"):
                        vuetify.VSelect(
                            label="Max Level",
                            v_model=("max_level",),
                            items=(
                                generalFunctions.get_default(
                                    "max_level_list", "default_values"
                                ),
                            ),
                            dense=True,
                        )
                with vuetify.VCol(classes="pa-0"):
                    vuetify.VListItemSubtitle(
                        "nCell",
                        classes="font-weight-bold black--text",
                    )
                with vuetify.VRow(classes="my-0"):
                    for direction in ["x", "y", "z"]:
                        with vuetify.VCol(cols=4, classes="py-0"):
                            vuetify.VTextField(
                                placeholder=direction,
                                v_model=(f"n_cell_{direction}",),
                                error_messages=(f"error_message_n_cell_{direction}",),
                                type="number",
                                step=generalFunctions.get_default("n_cell", "steps"),
                                __properties=["step"],
                                dense=True,
                                style="margin-top: -5px",
                            )
                with vuetify.VCol(classes="pa-0"):
                    vuetify.VListItemSubtitle(
                        "Blocking Factor",
                        classes="font-weight-bold black--text mt-2",
                    )
                with vuetify.VRow(classes="my-0"):
                    for direction in ["x", "y", "z"]:
                        with vuetify.VCol(cols=4, classes="py-0"):
                            vuetify.VTextField(
                                placeholder=direction,
                                v_model=(f"blocking_factor_{direction}",),
                                error_messages=(
                                    f"error_message_blocking_factor_{direction}",
                                ),
                                type="number",
                                step=generalFunctions.get_default(
                                    "blocking_factor", "steps"
                                ),
                                __properties=["step"],
                                dense=True,
                                style="margin-top: -5px",
                            )
                with vuetify.VCol(classes="pa-0"):
                    vuetify.VListItemSubtitle(
                        "prob_relative",
                        classes="font-weight-bold black--text mt-2",
                    )
                with vuetify.VRow(classes="my-0"):
                    with vuetify.VCol(
                        v_for=("(field, index) in prob_relative_fields",),
                        classes="py-0",
                    ):
                        vuetify.VTextField(
                            placeholder=("val."),
                            v_model=("field.value",),
                            input=(ctrl.update_prob_relative, "[index, $event]"),
                            error_messages=("field.error_message",),
                            type="number",
                            step=("field.step",),
                            __properties=["step"],
                            dense=True,
                            style="margin-top: -5px",
                        )

    @staticmethod
    def dialog_space_charge_settings():
        """
        Creates UI content for space charge configuration
        settings.
        """
        with vuetify.VCard():
            with vuetify.VTabs(
                v_model=("space_charge_tab", "Advanced Multigrid Settings")
            ):
                vuetify.VTab("Settings")
            vuetify.VDivider()
            with vuetify.VTabsItems(v_model="space_charge_tab"):
                with vuetify.VTabItem():
                    with vuetify.VContainer(fluid=True):
                        with vuetify.VRow(
                            classes="my-2", v_if="poisson_solver == 'multigrid'"
                        ):
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Relative Tolerance",
                                    v_model=("mlmg_relative_tolerance",),
                                    error_messages=(
                                        "error_message_mlmg_relative_tolerance",
                                    ),
                                    type="number",
                                    step=generalFunctions.get_default(
                                        "mlmg_relative_tolerance", "steps"
                                    ),
                                    __properties=["step"],
                                    dense=True,
                                )
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Absolute Tolerance",
                                    v_model=("mlmg_absolute_tolerance",),
                                    error_messages=(
                                        "error_message_mlmg_absolute_tolerance",
                                    ),
                                    suffix=generalFunctions.get_default(
                                        "mlmg_absolute_tolerance", "units"
                                    ),
                                    type="number",
                                    step=generalFunctions.get_default(
                                        "mlmg_absolute_tolerance", "steps"
                                    ),
                                    __properties=["step"],
                                    dense=True,
                                )
                        with vuetify.VRow(
                            classes="my-0", v_if="poisson_solver == 'multigrid'"
                        ):
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Max Iterations",
                                    v_model=("mlmg_max_iters",),
                                    error_messages=("error_message_mlmg_max_iters",),
                                    type="number",
                                    step=generalFunctions.get_default(
                                        "mlmg_max_iters", "steps"
                                    ),
                                    __properties=["step"],
                                    dense=True,
                                )
                            with vuetify.VCol(cols=6, classes="py-0"):
                                vuetify.VTextField(
                                    label="MLMG Verbosity",
                                    v_model=("mlmg_verbosity",),
                                    error_messages=("error_message_mlmg_verbosity",),
                                    type="number",
                                    step=generalFunctions.get_default(
                                        "mlmg_verbosity", "steps"
                                    ),
                                    __properties=["step"],
                                    dense=True,
                                )
