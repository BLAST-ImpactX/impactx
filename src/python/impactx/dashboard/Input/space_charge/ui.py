"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import ctrl, state, vuetify
from ...Input.components import (
    CardBase,
    CardComponents,
    InputComponents,
    NavigationComponents,
)
from .. import GeneralFunctions
from ..validation import DashboardValidation, errors_tracker
from .utils import SpaceChargeFunctions

state.prob_relative = []
state.prob_relative_fields = []
state.n_cell = []


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------


def populate_prob_relative_fields():
    tot_num_prob_relative_fields = int(state.max_level) + 1
    fft_first_field_value = GeneralFunctions.get_default(
        "prob_relative_first_value_fft", "default_values"
    )
    multigrid_first_field_value = GeneralFunctions.get_default(
        "prob_relative_first_value_multigrid", "default_values"
    )
    first_field_value = 0
    if state.poisson_solver == "fft":
        first_field_value = fft_first_field_value
    elif state.poisson_solver == "multigrid":
        first_field_value = multigrid_first_field_value

    if state.prob_relative:
        size = len(state.prob_relative)
        num_of_extra_fields = tot_num_prob_relative_fields - size

        if size < tot_num_prob_relative_fields:
            state.prob_relative.extend([0.0] * (num_of_extra_fields))
        elif size > tot_num_prob_relative_fields:
            state.prob_relative = state.prob_relative[:tot_num_prob_relative_fields]
    else:
        state.prob_relative = [first_field_value] + [0.0] * (
            tot_num_prob_relative_fields - 1
        )

    state.prob_relative_fields = [
        {
            "value": state.prob_relative[i],
            "error_message": DashboardValidation.validate_prob_relative_fields(
                i, state.prob_relative[i]
            ),
            "step": GeneralFunctions.get_default("prob_relative", "steps"),
        }
        for i in range(tot_num_prob_relative_fields)
    ]


# -----------------------------------------------------------------------------
# Decorators
# -----------------------------------------------------------------------------
@state.change("poisson_solver")
def on_poisson_solver_change(poisson_solver, **kwargs):
    populate_prob_relative_fields()
    state.dirty("prob_relative_fields")
    errors_tracker.update_simulation_validation_status()


@state.change("space_charge")
def on_space_charge_change(space_charge, **kwargs):
    state.dynamic_size = space_charge != "false"
    errors_tracker.update_simulation_validation_status()


@state.change("max_level")
def on_max_level_change(max_level, **kwargs):
    populate_prob_relative_fields()
    errors_tracker.update_simulation_validation_status()


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
    error_message = DashboardValidation.validate_prob_relative_fields(
        index, prob_relative_value
    )

    state.prob_relative_fields[index]["value"] = value
    state.prob_relative_fields[index]["error_message"] = error_message

    # Validate next index if it exists
    if index + 1 < len(state.prob_relative):
        next_value = state.prob_relative[index + 1]

        next_error_message = DashboardValidation.validate_prob_relative_fields(
            index + 1, next_value
        )
        state.prob_relative_fields[index + 1]["error_message"] = next_error_message

    state.dirty("prob_relative_fields")
    errors_tracker.update_simulation_validation_status()


# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------


class SpaceChargeConfiguration(CardBase):
    HEADER_NAME = "Space Charge"

    def __init__(self):
        super().__init__()

    def card_content(self):
        with vuetify.VDialog(
            v_model=("space_charge_dialog_settings", False), width="500px"
        ):
            SpaceChargeConfiguration.dialog_settings()

        with vuetify.VCard(**self.card_props):
            CardComponents.input_header(
                self.HEADER_NAME,
                additional_components={"end": SpaceChargeFunctions.multigrid_settings},
            )
            with vuetify.VCardText(**self.CARD_TEXT_OVERFLOW):
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(cols=4):
                        InputComponents.select(
                            label="Poisson Solver",
                        )
                    with vuetify.VCol(cols=4):
                        InputComponents.select(
                            label="Particle Shape",
                        )
                    with vuetify.VCol(cols=4):
                        InputComponents.select(
                            label="Max Level",
                        )
                for field in ["n_cell", "blocking_factor"]:
                    with vuetify.VRow(**self.ROW_STYLE):
                        for direction in ["x", "y", "z"]:
                            with vuetify.VCol(cols=4):
                                InputComponents.text_field(
                                    label=field if direction == "x" else "",
                                    v_model_name=f"{field}_{direction}",
                                    prefix=f"{direction}:",
                                )
                with vuetify.VRow(**self.ROW_STYLE):
                    with vuetify.VCol(
                        cols=4,
                        v_for="(field, index) in prob_relative_fields",
                        **self.ROW_STYLE,
                    ):
                        vuetify.VTextField(
                            label=("index === 0 ? 'prob_relative' : ''",),
                            v_model=("field.value",),
                            update_modelValue=(
                                ctrl.update_prob_relative,
                                "[index, $event]",
                            ),
                            error_messages=("field.error_message",),
                            type="number",
                            step=("field.step",),
                            __properties=["step"],
                            density="compact",
                            variant="underlined",
                            hide_details="auto",
                        )

    @staticmethod
    def dialog_settings():
        """
        Creates UI content for space charge configuration
        settings.
        """
        dialog_name = "space_charge_dialog_tab_settings"
        with NavigationComponents.create_dialog_tabs(
            dialog_name, 1, ["Advanced Multigrid Settings"]
        ):
            with vuetify.VTabsWindow(v_model=("dialog_name", 0)):
                with vuetify.VTabsWindowItem():
                    with vuetify.VCardText():
                        with vuetify.VRow():
                            with vuetify.VCol():
                                InputComponents.text_field(
                                    label="MLMG Relative Tolerance",
                                )
                            with vuetify.VCol():
                                InputComponents.text_field(
                                    label="MLMG Absolute Tolerance",
                                )
                        with vuetify.VRow():
                            with vuetify.VCol():
                                InputComponents.text_field(
                                    label="MLMG Max Iters",
                                )
                            with vuetify.VCol():
                                InputComponents.text_field(
                                    label="MLMG Verbosity",
                                )
