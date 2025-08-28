"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from .... import ctrl, html, state, vuetify
from . import components
from . import utils as _utils

_INIT_VALUE = ""
state.lattice_defaults = [
    {"name": _INIT_VALUE, "value": _INIT_VALUE, "error_message": _INIT_VALUE}
]
state.is_only_default = len(state.lattice_defaults) == 1
state.lattice_defaults_filter = ""
state.lattice_defaults_page = 1
state.lattice_defaults_filtered = []
state.lattice_defaults_no_results = False


# -----------------------------------------------------------------------------
# State listeners and controllers
# -----------------------------------------------------------------------------


@state.change("lattice_defaults")
def _on_defaults_change():
    _utils.update_delete_availability()
    _utils.apply_overrides_to_parameter_map()
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


@ctrl.add("reset_lattice_defaults")
def _on_reset_defaults() -> None:
    state.lattice_defaults = _utils.build_initial_defaults_list()
    state.lattice_defaults_filter = ""
    state.dirty("lattice_defaults")


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
        with vuetify.VCardText(classes="py-2"):
            components.search_bar()
        with vuetify.VCardText():
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
            components.pagination()
            with vuetify.VRow(classes="mt-2"):
                with vuetify.VCol():
                    components.reset_button()
