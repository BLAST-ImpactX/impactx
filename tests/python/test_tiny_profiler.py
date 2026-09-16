#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Tests for the default AMReX tiny profiler report location.

By default ImpactX writes the tiny profiler report to
``<diag_file_prefix>/performance.txt`` (e.g. ``diags/performance.txt``) instead of
stdout, and turns the profiler off when diagnostics are disabled so that no
diagnostics folder is created.

The shared pytest ``conftest.py`` disables the tiny profiler process-wide
(``tiny_profiler.enabled=0``), so these tests manage their own AMReX lifecycle
(``manages_amrex``) and enable the profiler explicitly. AMReX re-reads
``tiny_profiler.output_file`` once per init/finalize cycle, so each test resolves
its own output location independently.
"""

from pathlib import Path

import pytest

from impactx import ImpactX, distribution, elements


def run_minimal_simulation(sim):
    """Run a minimal drift simulation on an already-initialized sim."""
    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    distr = distribution.Waterbag(
        lambdaX=1.0e-3,
        lambdaY=1.0e-3,
        lambdaT=1.0e-3,
        lambdaPx=1.0e-4,
        lambdaPy=1.0e-4,
        lambdaPt=1.0e-3,
    )
    sim.add_particles(1.0e-9, distr, 100)

    sim.lattice.append(elements.Drift(name="d1", ds=0.25))
    sim.track_particles()


@pytest.mark.manages_amrex
def test_performance_file_default():
    """By default, the profiler report is written to diags/performance.txt."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.tiny_profiler = True
    sim.init_grids()

    run_minimal_simulation(sim)

    sim.finalize()

    performance = Path("diags") / "performance.txt"
    assert performance.exists()
    assert performance.stat().st_size > 0


@pytest.mark.manages_amrex
def test_performance_file_custom_prefix():
    """The profiler report follows a custom diag_file_prefix."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.tiny_profiler = True
    sim.diag_file_prefix = "perf_out"
    sim.init_grids()

    run_minimal_simulation(sim)

    sim.finalize()

    assert (Path("perf_out") / "performance.txt").exists()
    assert not Path("diags").exists()


@pytest.mark.manages_amrex
def test_performance_file_diagnostics_disabled():
    """With diagnostics disabled, the profiler is off and no diags folder is created."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.tiny_profiler = True
    sim.diagnostics = False
    sim.init_grids()

    run_minimal_simulation(sim)

    sim.finalize()

    assert not Path("diags").exists()
    assert not Path("performance.txt").exists()


@pytest.mark.manages_amrex
def test_performance_file_explicit_override():
    """An explicit tiny_profiler_file takes precedence over the default."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.tiny_profiler = True
    sim.tiny_profiler_file = "myperf.txt"
    sim.init_grids()

    run_minimal_simulation(sim)

    sim.finalize()

    assert Path("myperf.txt").exists()
    assert not (Path("diags") / "performance.txt").exists()
