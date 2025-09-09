"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from .... import ctrl, html, state, vuetify
from copy import deepcopy
from . import components
from . import utils as _utils

# Initialize applied and staged defaults
state.lattice_defaults_applied = _utils.build_initial_defaults_list()
state.lattice_defaults = deepcopy(state.lattice_defaults_applied)
state.is_only_default = len(state.lattice_defaults) == 1
state.lattice_defaults_filter = ""
state.lattice_defaults_page = 1
state.lattice_defaults_filtered = []
state.lattice_defaults_no_results = False
state.lattice_defaults_has_changes = False


def _sync_has_changes_flag() -> None:
    """
    Set has_changes when staged edits differ from applied values.
    """
    state.lattice_defaults_has_changes = (
        state.lattice_defaults != state.lattice_defaults_applied
    )
    state.dirty("lattice_defaults_has_changes")


# -----------------------------------------------------------------------------
# State listeners and controllers
# -----------------------------------------------------------------------------


@state.change("lattice_defaults")
def _on_defaults_change(*_args, **_kwargs):
    # Only update UI-related state; do not apply globally until "Apply" is clicked
    _utils.update_delete_availability()
    _utils.sync_filtered_defaults()


@state.change("lattice_defaults_filter")
def _on_filter_change(*_args, **_kwargs):
    """
    Reset to the first page when filter changes and update filtered list
    """
    state.lattice_defaults_page = 1
    state.dirty("lattice_defaults_page")
    _utils.sync_filtered_defaults()


@ctrl.add("update_lattice_default")
def _on_update_default(field_name: str, identifier, new_value) -> None:
    """
    Update a field of a lattice default entry.
    field_name: which field to update, e.g. "value"
    identifier: either the row index (int) or the parameter name (str)
    new_value: the value to assign
    """
    # Resolve the index: use it directly if it's int, otherwise look up by name
    if isinstance(identifier, int):
        row_index = identifier
    else:
        row_index = next(
            (
                i
                for i, row in enumerate(state.lattice_defaults)
                if row.get("name") == identifier
            ),
            None,
        )
        if row_index is None:
            return

    entry = state.lattice_defaults[row_index]
    if field_name == "value":
        entry["value"] = new_value

    state.dirty("lattice_defaults")
    _sync_has_changes_flag()


@ctrl.add("reset_lattice_defaults")
def _on_reset_defaults() -> None:
    state.lattice_defaults = _utils.build_initial_defaults_list()
    state.lattice_defaults_filter = ""
    state.dirty("lattice_defaults")
    _sync_has_changes_flag()


@ctrl.add("apply_lattice_defaults")
def _on_apply_defaults() -> None:
    """
    Mark current overrides as applied. Recompute parameter map to ensure
    downstream consumers see the latest values, then clear the dirty flag.
    """
    # Persist staged edits as applied and recompute parameter map
    state.lattice_defaults_applied = deepcopy(state.lattice_defaults)
    _utils.apply_overrides_to_parameter_map()
    state.lattice_defaults_has_changes = False
    state.dirty("lattice_defaults_has_changes")


@state.change("lattice_configuration_dialog_settings")
def _on_dialog_visibility_change(lattice_configuration_dialog_settings, **_):
    # On close/cancel, revert staged edits to last applied and clear dirty flag
    if not lattice_configuration_dialog_settings:
        state.lattice_defaults = deepcopy(state.lattice_defaults_applied)
        _utils.sync_filtered_defaults()
        state.lattice_defaults_has_changes = False
        state.dirty("lattice_defaults", "lattice_defaults_has_changes")


class LatticeDefaultsHandler:
    """
    UI entry point for editing lattice parameter default overrides.
    """

    @staticmethod
    def _build_initial_defaults_list():
        return _utils.build_initial_defaults_list()

    @staticmethod
    def defaults_handler():
        """
        Renders the Defaults tab using a Variables-like UI.
        """
        if (
            isinstance(state.lattice_defaults, list)
            and len(state.lattice_defaults) == 1
            and not (
                state.lattice_defaults[0].get("name")
                or state.lattice_defaults[0].get("value")
            )
        ):
            state.lattice_defaults = _utils.build_initial_defaults_list()
            state.dirty("lattice_defaults")
            _utils.sync_filtered_defaults()
        elif not state.lattice_defaults_filtered:
            _utils.sync_filtered_defaults()
        with vuetify.VCardText(classes="py-1 pb-0"):
            components.search_bar()
        with vuetify.VCardText(classes="pt-1"):
            with vuetify.VContainer(fluid=True):
                with vuetify.VRow(
                    v_for=(
                        "(item, index) in lattice_defaults_filtered.slice((lattice_defaults_page - 1) * 5, (lattice_defaults_page) * 5)",
                    ),
                    classes="align-center justify-center py-0",
                ):
                    with vuetify.VCol(cols=5, classes="pr-0"):
                        components.text_field(
                            placeholder="Parameter Name",
                            v_model=("item.name",),
                            id=("'default_name_' + (item.name || '')",),
                            background_color="grey lighten-4",
                            readonly=True,
                            clearable=False,
                        )
                    with vuetify.VCol(cols=1, classes="px-0 text-center"):
                        html.Span("=", classes="mx-0")
                    with vuetify.VCol(cols=4, classes="pl-0"):
                        components.text_field(
                            placeholder="Default Value",
                            v_model=("item.value",),
                            id=("'default_value_' + (item.name || '')",),
                            type="text",
                            background_color="grey lighten-4",
                            update_modelValue=(
                                ctrl.update_lattice_default,
                                "['value', item.name, $event]",
                            ),
                            clearable=True,
                        )
                vuetify.VAlert(
                    "No matches found",
                    type="info",
                    variant="tonal",
                    density="compact",
                    border=True,
                    classes="ma-2",
                    v_show=("lattice_defaults_no_results",),
                )
        with vuetify.VCardText(classes="pt-0 pb-2"):
            # Pagination centered above the footer actions
            with vuetify.VRow(classes="justify-center"):
                with vuetify.VCol(cols="auto"):
                    components.pagination()

            vuetify.VDivider(classes="my-2")

            # Footer actions: Reset (left) and Close/Apply (right)
            with vuetify.VRow(classes="align-center"):
                with vuetify.VCol(cols=6):
                    vuetify.VBtn(
                        "Reset",
                        color="primary",
                        variant="tonal",
                        click=ctrl.reset_lattice_defaults,
                    )
                with vuetify.VCol(cols=6, classes="d-flex justify-end"):
                    vuetify.VBtn(
                        "Close",
                        variant="text",
                        color="#00313C",
                        click="lattice_configuration_dialog_settings = false",
                        classes="mr-2",
                    )
                    vuetify.VBtn(
                        "Apply",
                        color="primary",
                        disabled=("!lattice_defaults_has_changes",),
                        click=ctrl.apply_lattice_defaults,
                    )
