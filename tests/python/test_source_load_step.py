#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Test that the Source element can load a selected step (openPMD iteration).
"""

import re
from pathlib import Path

import pytest

from impactx import Config, ImpactX, distribution, elements, push

io = pytest.importorskip("openpmd_api")

if not Config.have_openpmd:
    pytest.skip("ImpactX was compiled without openPMD support", allow_module_level=True)

# the drift between the two beam monitors: long enough that the beam widens
# well beyond the sampling noise of the moments of a finite number of particles
DRIFT_DS = 2.0

# HDF5 as in the solenoid_restart example, else whatever this build provides
BACKEND = "h5" if Config.openpmd_backends.get("hdf5", False) else "default"

# exact for a single drift, but single precision accumulates a few epsilon
RTOL = 1.0e-12 if Config.precision == "DOUBLE" else 1.0e-5


def write_series(name, npart, periods=1):
    """Track a beam through two monitors and return the written series' path.

    The lattice writes two steps per period: one at s=0 and one at s=DRIFT_DS.
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.periods = periods
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(250.0)

    distr = distribution.Waterbag(
        lambdaX=1.559531175539e-3,
        lambdaY=2.205510139392e-3,
        lambdaT=1.0e-3,
        lambdaPx=6.41218345413e-4,
        lambdaPy=9.06819680526e-4,
        lambdaPt=1.0e-3,
    )
    sim.add_particles(1.0e-9, distr, npart)

    monitor = elements.BeamMonitor(name, backend=BACKEND)
    sim.lattice.extend([monitor, elements.Drift(name="d1", ds=DRIFT_DS), monitor])

    try:
        sim.track_particles()
    finally:
        # this closes the openPMD series, so that we can read it back below
        sim.finalize()

    files = sorted(Path("diags/openPMD").glob(f"{name}.*"))
    assert files, f"no openPMD series found for monitor '{name}'"

    return str(files[0])


def stored_steps(series_path):
    """The steps (openPMD iterations) stored in a series."""
    series = io.Series(series_path, io.Access.read_only)
    steps = sorted(series.iterations)
    series.close()

    return steps


@pytest.mark.manages_amrex
def test_source_load_step():
    """
    This tests that the Source element loads the step (openPMD iteration)
    selected by step number with load_step and by position with load_step_index.
    """
    npart = 512
    series_path = write_series("mon_load_step", npart)

    steps = stored_steps(series_path)
    assert len(steps) == 2

    # read back into a fresh simulation
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    # keep init_grids from moving the diags directory we just wrote out of the way
    sim.diagnostics = False
    sim.init_grids()
    beam = sim.beam

    try:
        # sig_x of the loaded particles, per step: the reference particle alone
        # would not show whether the particles come from the selected step
        sig_x = {}

        for selection, s_expected in [
            ({"load_step": steps[0]}, 0.0),  # first step: by step number
            ({"load_step": steps[-1]}, DRIFT_DS),  # last step: by step number
            ({"load_step_index": 0}, 0.0),  # first step: by position
            ({"load_step_index": 1}, DRIFT_DS),  # last step: by position
            ({"load_step_index": -1}, DRIFT_DS),  # last step: counted back
            ({"load_step_index": -2}, 0.0),  # first step: counted back
            ({}, DRIFT_DS),  # the default is the last step in the file
        ]:
            # the source adds to the container, so read each step into an empty one
            beam.clear_particles()

            source = elements.Source("openPMD", series_path, name="source", **selection)
            push(beam, source)
            assert beam.ref.s == pytest.approx(s_expected, rel=RTOL, abs=1.0e-12)
            assert beam.ref.kin_energy_MeV == pytest.approx(250.0, rel=RTOL)

            moments = beam.beam_moments()
            # the whole beam was read into the emptied container
            assert moments["charge_C"] == pytest.approx(1.0e-9, rel=RTOL)
            sig_x.setdefault(s_expected, []).append(moments["sig_x"])

        # the same step selected by number or by position holds the same particles
        for values in sig_x.values():
            assert values == pytest.approx([values[0]] * len(values), rel=RTOL)

        # ... and the beam widens along the drift between the two steps
        assert sig_x[DRIFT_DS][0] > 1.1 * sig_x[0.0][0]

        # a step that is not in the file lists the steps that are in it
        missing = steps[-1] + 1
        with pytest.raises(RuntimeError, match="available steps") as error:
            push(beam, elements.Source("openPMD", series_path, load_step=missing))
        assert f"{len(steps)} available steps" in str(error.value)

        # an index that reaches past either end of the file
        for load_step_index in [-3, 2]:
            with pytest.raises(RuntimeError, match="out of range") as error:
                push(
                    beam,
                    elements.Source(
                        "openPMD", series_path, load_step_index=load_step_index
                    ),
                )
            assert f"{len(steps)} available steps" in str(error.value)

        # both options set after construction is caught when the step is read
        source = elements.Source("openPMD", series_path, load_step=steps[0])
        source.load_step_index = -1
        with pytest.raises(ValueError, match="not both"):
            push(beam, source)
    finally:
        sim.finalize()


def test_source_load_step_selection():
    """
    This tests that at most one of load_step and load_step_index selects the
    step and that load_step is a step number, not a position in the file.
    """
    with pytest.raises(ValueError, match="not both"):
        elements.Source("openPMD", "beam.h5", load_step=3, load_step_index=-1)

    # counting back from the last step is load_step_index, not a negative load_step
    with pytest.raises(ValueError, match="load_step_index"):
        elements.Source("openPMD", "beam.h5", load_step=-2)


@pytest.mark.manages_amrex
def test_source_load_step_many_steps():
    """
    This tests that the available steps in the error message are elided in the
    middle for a series with many steps, e.g. one per turn in a ring.
    """
    series_path = write_series("mon_load_step_many", npart=32, periods=51)

    steps = stored_steps(series_path)
    assert len(steps) > 100

    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    # keep init_grids from moving the diags directory we just wrote out of the way
    sim.diagnostics = False
    sim.init_grids()
    beam = sim.beam

    try:
        missing = steps[-1] + 1
        with pytest.raises(RuntimeError, match="available steps") as error:
            push(beam, elements.Source("openPMD", series_path, load_step=missing))

        message = str(error.value)
        assert "..." in message
        assert f"{len(steps)} available steps" in message

        # the first and last steps are listed, the middle ones are elided
        listed = re.findall(r"\d+", message.split("available steps: ")[1])
        assert str(steps[0]) in listed
        assert str(steps[-1]) in listed
        assert str(steps[len(steps) // 2]) not in listed
        assert len(listed) == 100
    finally:
        sim.finalize()


def test_source_load_step_serialization():
    """
    This tests that load_step and load_step_index are part of the element's
    repr and to_dict(), and that an unset option is None.
    """
    source = elements.Source("openPMD", "beam.h5", load_step=3)
    assert source.load_step == 3
    assert source.load_step_index is None
    assert "load_step=3" in repr(source)
    assert "load_step_index=None" in repr(source)

    d = source.to_dict()
    assert d["load_step"] == 3
    assert d["load_step_index"] is None

    # the dict round-trips through the constructor
    kwargs = {k: v for k, v in d.items() if k not in ("type", "ds")}
    assert elements.Source(**kwargs).load_step == 3

    source = elements.Source("openPMD", "beam.h5", load_step_index=-2)
    assert source.load_step is None
    assert source.load_step_index == -2
    assert source.to_dict()["load_step_index"] == -2

    # neither is set by default: the last step in the file is loaded
    default = elements.Source("openPMD", "beam.h5")
    assert default.load_step is None
    assert default.load_step_index is None
    assert "load_step=None" in repr(default)

    # an option can be unset again
    source.load_step_index = None
    assert source.load_step_index is None
