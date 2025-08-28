"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from impactx import elements

from .... import ctrl, state
from ...defaults_helper import InputDefaultsHelper
from ...utils import GeneralFunctions


EXCLUDED_LATTICE_DEFAULTS: set[str] = {"name"}
EXCLUDED_LATTICE_CLASSES: set[str] = {"BeamMonitor"}

def apply_overrides_to_parameter_map() -> None:
    """
    Rebuild the element->(name, default, type) map and apply overrides.
    Cast numeric overrides to the appropriate type for downstream use.
    """
    data = InputDefaultsHelper.class_parameters_with_defaults(elements)

    # Quick lookup of parameter types for casting
    type_lookup: dict[str, str] = {}
    for params in data.values():
        for pname, _pdefault, ptype in params:
            type_lookup[pname] = ptype

    for override in state.lattice_defaults:
        name = (override.get("name") or "").strip()
        value = override.get("value")
        if not name:
            continue
        if name in EXCLUDED_LATTICE_DEFAULTS:
            continue

        ptype = type_lookup.get(name)
        cast_value = value
        if ptype in {"int", "float"}:
            cast_value = GeneralFunctions.convert_to_numeric(value)

        for key, parameters in data.items():
            if key in EXCLUDED_LATTICE_CLASSES:
                continue
            for i, (pname, _pdefault, ptype) in enumerate(parameters):
                if pname == name:
                    parameters[i] = (pname, cast_value, ptype)

    state.listOfLatticeElementParametersAndDefault = data


def update_delete_availability() -> None:
    state.is_only_default = len(state.lattice_defaults) == 1


def build_initial_defaults_list() -> list[dict]:
    """
    Build an initial list of defaults exposing all known lattice parameter names
    with their current default values.
    """
    rows: list[dict] = []
    seen: set[str] = set()
    data = InputDefaultsHelper.class_parameters_with_defaults(elements)
    for class_name, params in data.items():
        if class_name in EXCLUDED_LATTICE_CLASSES:
            continue
        for pname, pdefault, _ptype in params:
            if pname in EXCLUDED_LATTICE_DEFAULTS:
                continue
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


def _filter_defaults_by_query() -> list[dict]:
    query = (state.lattice_defaults_filter or "").lower()
    defaults = state.lattice_defaults or []
    if not query:
        return defaults
    return [row for row in defaults if query in (row.get("name") or "").lower()]


def sync_filtered_defaults() -> None:
    filtered = _filter_defaults_by_query()
    state.lattice_defaults_filtered = filtered
    state.lattice_defaults_no_results = bool(state.lattice_defaults_filter) and not filtered
    state.dirty("lattice_defaults_filtered", "lattice_defaults_no_results")