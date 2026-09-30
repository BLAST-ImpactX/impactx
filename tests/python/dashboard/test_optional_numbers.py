"""
This file is part of ImpactX

Copyright 2022-2026 The ImpactX Community
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

from impactx import elements
from impactx.dashboard import state
from impactx.dashboard.Input.defaults_helper import InputDefaultsHelper
from impactx.dashboard.Input.validation import DashboardValidation
from impactx.dashboard.Input.validation.inputs import FLOAT_ERROR_MESSAGE
from impactx.dashboard.Run.simulation import build_lattice_list
from impactx.dashboard.Toolbar.file_imports import ui_populator


def test_optional_numbers_are_parsed_from_the_signature():
    """A bend's rc, phi and B default to None: optional floats"""
    parameters = InputDefaultsHelper.extract_parameters(elements.Sbend.__init__.__doc__)
    optional = {name: (default, kind) for name, default, kind in parameters}
    for name in ["rc", "phi", "B"]:
        assert optional[name] == ("None", "float | None"), name
    assert optional["ds"] == (None, "float")


def test_optional_numbers_may_be_left_unset():
    for unset in [None, "None", "", "  "]:
        assert (
            DashboardValidation.validate(
                "rc", unset, category="lattice", parameter_type="float | None"
            )
            == []
        )
    assert (
        DashboardValidation.validate(
            "rc", "10.5", category="lattice", parameter_type="float | None"
        )
        == []
    )
    assert DashboardValidation.validate(
        "rc", "ten", category="lattice", parameter_type="float | None"
    ) == [FLOAT_ERROR_MESSAGE]


def test_unset_optional_numbers_are_exported_as_none():
    previous = state.selected_lattice_list
    state.selected_lattice_list = [
        {
            "name": "Sbend",
            "parameters": [
                {"parameter_name": "ds", "sim_input": "0.5", "parameter_type": "float"},
                {
                    "parameter_name": "rc",
                    "sim_input": "10.0",
                    "parameter_type": "float | None",
                },
                {
                    "parameter_name": "phi",
                    "sim_input": "None",
                    "parameter_type": "float | None",
                },
                {
                    "parameter_name": "B",
                    "sim_input": "",
                    "parameter_type": "float | None",
                },
            ],
        }
    ]
    try:
        assert (
            "elements.Sbend(ds=0.5, rc=10.0, phi=None, B=None)" in build_lattice_list()
        )
    finally:
        state.selected_lattice_list = previous


def test_a_file_is_imported_once(monkeypatch):
    """A state flush during the import re-runs the import listener; it must not re-import."""
    imports = []

    def populate(file):
        imports.append(file)
        # what trame does when the import flushes the state: run the listener again
        ui_populator.on_import_file_change(file)

    monkeypatch.setattr(ui_populator.DashboardParser, "file_details", lambda file: None)
    monkeypatch.setattr(
        ui_populator, "populate_impactx_simulation_file_to_ui", populate
    )

    ui_populator.on_import_file_change({"name": "run.py", "content": b""})
    assert len(imports) == 1

    # the next file is imported again
    ui_populator.on_import_file_change({"name": "run2.py", "content": b""})
    assert len(imports) == 2
