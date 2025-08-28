from typing import Optional, Tuple

from impactx import elements

from ... import ctrl, vuetify, html, state
from ..defaults_helper import InputDefaultsHelper
from ..utils import GeneralFunctions


init_value = ""
state.lattice_defaults = [
    {"name": init_value, "value": init_value, "error_message": init_value}
]
state.is_only_default = len(state.lattice_defaults) == 1
state.lattice_defaults_filter = ""
state.lattice_defaults_page = 1
# fixed page size of 5 rows per page
state.lattice_defaults_filtered = []
state.lattice_defaults_no_results = False


class LatticeDefaultsHandler:
    """
    UI and logic for editing lattice parameter default overrides.

    - Presents a list of (name, value) rows, similar to the variables UI.
    - Validates names against known lattice parameter names.
    - Applies overrides to `state.listOfLatticeElementParametersAndDefault` so new
      lattice elements use these values as the initial defaults.
    """

    # -------------------------------------------------------------------------
    # Internal helpers
    # -------------------------------------------------------------------------
    @staticmethod
    def _known_parameter_names() -> set[str]:
        names = set()
        data = InputDefaultsHelper.class_parameters_with_defaults(elements)
        for params in data.values():
            for pname, _pdefault, _ptype in params:
                names.add(pname)
        return names

    @staticmethod
    def _apply_overrides_to_parameter_map() -> None:
        """
        Rebuilds the element->(name, default, type) map and applies overrides.
        Casts numeric overrides to the appropriate type.
        """
        data = InputDefaultsHelper.class_parameters_with_defaults(elements)

        # Make a quick lookup of parameter types for casting
        type_lookup: dict[str, str] = {}
        for params in data.values():
            for pname, _pdefault, ptype in params:
                type_lookup[pname] = ptype

        for override in state.lattice_defaults:
            name = (override.get("name") or "").strip()
            value = override.get("value")
            if not name:
                continue

            ptype = type_lookup.get(name)
            cast_value = value
            if ptype in {"int", "float"}:
                cast_value = GeneralFunctions.convert_to_numeric(value)

            for key, parameters in data.items():
                for i, (pname, _pdefault, ptype) in enumerate(parameters):
                    if pname == name:
                        parameters[i] = (pname, cast_value, ptype)

        state.listOfLatticeElementParametersAndDefault = data

    @staticmethod
    def _update_delete_availability() -> None:
        state.is_only_default = len(state.lattice_defaults) == 1

    @staticmethod
    def _duplicate_indexes(new_name: str, current_index: int) -> list[int]:
        dups = [
            idx
            for idx, item in enumerate(state.lattice_defaults)
            if item["name"] == new_name and idx != current_index
        ]
        if dups:
            dups.append(current_index)
        return dups

    @staticmethod
    def _build_initial_defaults_list() -> list[dict]:
        """
        Build an initial list of defaults exposing all known lattice parameter names
        with their current default values.
        """
        rows: list[dict] = []
        seen: set[str] = set()
        data = InputDefaultsHelper.class_parameters_with_defaults(elements)
        for params in data.values():
            for pname, pdefault, _ptype in params:
                if pname in seen:
                    continue
                seen.add(pname)
                rows.append(
                    {
                        "name": pname,
                        "value": pdefault if pdefault is not None else "",
                        "error_message": "",
                    }
                )
        return rows

    @staticmethod
    def _filter_defaults_by_query() -> list[dict]:
        """
        Returns defaults filtered by the current search query.
        """
        query = (state.lattice_defaults_filter or "").lower()
        defaults = state.lattice_defaults or []
        if not query:
            return defaults
        return [row for row in defaults if query in (row.get("name") or "").lower()]

    @staticmethod
    def _sync_filtered_defaults() -> None:
        filtered = LatticeDefaultsHandler._filter_defaults_by_query()
        state.lattice_defaults_filtered = filtered
        state.lattice_defaults_no_results = bool(state.lattice_defaults_filter) and not filtered
        state.dirty("lattice_defaults_filtered", "lattice_defaults_no_results")

    # -------------------------------------------------------------------------
    # State listeners
    # -------------------------------------------------------------------------
    @staticmethod
    @state.change("lattice_defaults")
    def on_defaults_change(lattice_defaults, **kwargs):
        LatticeDefaultsHandler._update_delete_availability()
        LatticeDefaultsHandler._apply_overrides_to_parameter_map()
        LatticeDefaultsHandler._sync_filtered_defaults()

    @staticmethod
    @state.change("lattice_defaults_filter")
    def on_filter_change(*_args, **_kwargs):
        # Reset to first page when filter changes and update filtered list
        state.lattice_defaults_page = 1
        state.dirty("lattice_defaults_page")
        LatticeDefaultsHandler._sync_filtered_defaults()

    # -------------------------------------------------------------------------
    # Controllers
    # -------------------------------------------------------------------------
    @staticmethod
    @ctrl.add("add_lattice_default")
    def on_add_default() -> None:
        # No-op: defaults list is fixed; names are read-only
        pass

    @staticmethod
    @ctrl.add("delete_lattice_default")
    def on_delete_default(index: int) -> None:
        # No-op: defaults list is fixed; names are read-only
        pass

    @staticmethod
    @ctrl.add("update_lattice_default")
    def on_update_default(key_name: str, key, event) -> None:
        # Resolve row index either directly or via parameter name
        idx = key
        if isinstance(key, str):
            idx = next(
                (
                    i
                    for i, row in enumerate(state.lattice_defaults)
                    if row.get("name") == key
                ),
                None,
            )
            if idx is None:
                return
        item = state.lattice_defaults[idx]
        if key_name == "value":
            item["value"] = event
        state.dirty("lattice_defaults")

    @staticmethod
    @ctrl.add("reset_lattice_defaults")
    def on_reset_defaults() -> None:
        state.lattice_defaults = LatticeDefaultsHandler._build_initial_defaults_list()
        state.lattice_defaults_filter = ""
        state.dirty("lattice_defaults")

    # -------------------------------------------------------------------------
    # UI
    # -------------------------------------------------------------------------
    @staticmethod
    def defaults_handler():
        """
        Renders the Defaults tab using a Variables-like UI.
        """
        # On first open, populate rows with all known defaults
        if (
            isinstance(state.lattice_defaults, list)
            and len(state.lattice_defaults) == 1
            and not (state.lattice_defaults[0].get("name") or state.lattice_defaults[0].get("value"))
        ):
            state.lattice_defaults = LatticeDefaultsHandler._build_initial_defaults_list()
            state.dirty("lattice_defaults")
            LatticeDefaultsHandler._sync_filtered_defaults()
        elif not state.lattice_defaults_filtered:
            # Ensure filtered list is initialized
            LatticeDefaultsHandler._sync_filtered_defaults()
        # Search + results summary (single line)
        with vuetify.VCardText(classes="py-2"):
            with vuetify.VRow(classes="align-center"):
                with vuetify.VCol(cols=True):
                    vuetify.VTextField(
                        v_model=("lattice_defaults_filter", ""),
                        label="Search parameters",
                        placeholder="e.g., nslice",
                        prepend_inner_icon="mdi-magnify",
                        variant="outlined",
                        density="compact",
                        hide_details=True,
                        clearable=True,
                        classes="text-body-2",
                        id="lattice_defaults_search",
                        aria_label="Search parameters",
                        __properties=["aria-label"],
                    )
                with vuetify.VCol(cols="auto", classes="d-flex align-center ml-3"):
                    vuetify.VIcon(
                        "mdi-checkbox-blank-circle",
                        size="x-small",
                        color="primary",
                        classes="mr-2",
                    )
                    html.Span(
                        "Matches: {{ lattice_defaults_filtered.length }} / {{ lattice_defaults.length }}",
                        classes="text-caption text-grey-darken-1",
                        aria_live="polite",
                        __properties=["aria-live"],
                    )

        # List area (no scrollbar, 5 items per page)
        with vuetify.VCardText():
            with vuetify.VContainer(fluid=True):
                with vuetify.VRow(
                    v_for=(
                        "(item, index) in lattice_defaults_filtered.slice((lattice_defaults_page - 1) * 5, (lattice_defaults_page) * 5)",
                    ),
                    classes="align-center justify-center py-0",
                ):
                    with vuetify.VCol(cols=5, classes="pr-0"):
                        vuetify.VTextField(
                            placeholder="Parameter Name",
                            v_model=("item.name",),
                            id=("'default_name_' + (item.name || '')",),
                            variant="outlined",
                            density="compact",
                            background_color="grey lighten-4",
                            readonly=True,
                            hide_details=True,
                            clearable=False,
                            style="min-width: 0;",
                        )
                    with vuetify.VCol(cols=1, classes="px-0 text-center"):
                        html.Span("=", classes="mx-0")
                    with vuetify.VCol(cols=4, classes="pl-0"):
                        vuetify.VTextField(
                            placeholder="Default Value",
                            v_model=("item.value",),
                            id=("'default_value_' + (item.name || '')",),
                            variant="outlined",
                            density="compact",
                            type="text",
                            background_color="grey lighten-4",
                            update_modelValue=(
                                ctrl.update_lattice_default,
                                "['value', item.name, $event]",
                            ),
                            hide_details=True,
                            clearable=True,
                            style="min-width: 0;",
                        )
                    with vuetify.VCol(cols=2, classes="d-flex"):
                        pass
                # No matches indicator
                vuetify.VAlert(
                    "No matches found",
                    type="info",
                    variant="tonal",
                    density="compact",
                    border=True,
                    classes="ma-2",
                    v_show=("lattice_defaults_no_results",),
                )

        # Pagination controls
        with vuetify.VCardText(classes="pt-0 pb-2"):
            vuetify.VPagination(
                v_model=("lattice_defaults_page", 1),
                length=("Math.max(1, Math.ceil(lattice_defaults_filtered.length / 5))",),
                total_visible=7,
                __properties=["length", "total_visible"],
                density="comfortable",
            )
            with vuetify.VRow(classes="mt-2"):
                with vuetify.VCol():
                    vuetify.VBtn(
                        "Reset Defaults",
                        id="reset_lattice_defaults",
                        color="primary",
                        click=ctrl.reset_lattice_defaults,
                        block=True,
                    )
