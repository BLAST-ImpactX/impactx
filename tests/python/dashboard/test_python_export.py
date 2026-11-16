"""
This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""

import time

from .utils import TIMEOUT


def exported_script(dashboard) -> str:
    """
    Returns the Python script that the Export button downloads.
    """
    return dashboard.sb.driver.execute_async_script(
        "const done = arguments[arguments.length - 1];"
        "window.trame.trigger('export').then(done, (error) => done(String(error)));"
    )


def assert_exported(dashboard, expected_lines: list) -> None:
    """
    Waits until the exported script contains all expected lines.
    """
    script = ""
    for _ in range(TIMEOUT):
        script = exported_script(dashboard)
        if all(line in script for line in expected_lines):
            return
        time.sleep(1)
    raise AssertionError(f"Exported script lacks {expected_lines}:\n{script}")


def test_python_export(dashboard):
    """
    The exported script sets what CSR and envelope tracking need.
    """
    # CSR without space charge needs the particle shape
    dashboard.set_input("csr", True)
    assert_exported(dashboard, ["sim.csr = True", "sim.particle_shape = 2"])

    # envelope tracking needs the beam intensity
    dashboard.set_input("tracking_mode", "Envelope Tracking")
    dashboard.set_input("space_charge", "3D")
    assert_exported(
        dashboard,
        ["bunch_charge_C = 1e-09", "sim.init_envelope(ref, distr, bunch_charge_C)"],
    )

    # 2D space charge takes the beam current instead of the bunch charge
    dashboard.set_input("space_charge", "2D")
    assert_exported(
        dashboard,
        ["beam_current_A = 0.5", "sim.init_envelope(ref, distr, beam_current_A)"],
    )
