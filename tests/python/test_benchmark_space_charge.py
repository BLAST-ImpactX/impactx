#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# This set of tests is used for performance benchmarking of the space-charge
# solvers (as micro-benchmarks). We use this file to rapidly evaluate
# performance changes when tuning the deposit -> Poisson-solve -> field-calc ->
# gather+push pipeline (``HandleSpacecharge``) on CPUs and GPUs.
#
# Unlike the per-element benchmarks in ``test_benchmark_elements.py``, the
# space-charge solve is not exposed as a standalone Python callable. We
# therefore drive the solve by tracking a short constant-focusing cell (a few
# slices) with space charge enabled, and benchmark ``sim.track_particles()``.
#
# -*- coding: utf-8 -*-

import os
from functools import partial

import pytest

from impactx import Config, ImpactX, distribution, elements

# benchmark config
if os.environ.get("IS_CODESPEED_CPU_SIMULATION") == "1":
    # https://codspeed.io/docs/instruments/cpu/index
    rounds = 1
    npart = 10_000
else:
    rounds = 5
    # space-charge solves are much heavier than a single element push (a full
    # Poisson solve over the mesh per round), so we use fewer particles than the
    # element benchmark to keep the (non-CodSpeed) ``ctest`` run reasonable.
    # Increase this for serious local benchmarking on capable hardware.
    npart = 100_000

# 2 GeV proton beam (matches the cfchannel example)
kin_energy_MeV = 2.0e3
bunch_charge_C = 1.0e-8


# Distinct space-charge solver code paths to benchmark. Each entry carries the
# configuration that must be set before ``init_grids()``. The ``fft`` flag marks
# cases that require ImpactX to be built with FFT support (IGF / 2D / 2p5D); on
# builds without it, those cases skip instead of erroring.
SC_CASES = [
    dict(
        id="3D_mlmg",
        space_charge="3D",
        poisson_solver="multigrid",
        n_cell=[32, 32, 32],
        blocking_factor=[16, 16, 16],
        prob_relative=[3.0],
        fft=False,
    ),
    dict(
        id="3D_fft",
        space_charge="3D",
        poisson_solver="fft",
        n_cell=[32, 32, 32],
        blocking_factor=[16, 16, 16],
        prob_relative=[1.1],
        fft=True,
    ),
    dict(
        id="2D_fft",
        space_charge="2D",
        poisson_solver="fft",
        n_cell=[256, 256, 1],
        blocking_factor=[128, 128, 1],
        prob_relative=[1.1],
        fft=True,
    ),
    dict(
        id="2p5D_fft",
        space_charge="2p5D",
        poisson_solver="fft",
        n_cell=[256, 256, 1],
        blocking_factor=[128, 128, 1],
        prob_relative=[1.1],
        fft=True,
    ),
    dict(
        id="Gauss3D",
        space_charge="Gauss3D",
        n_cell=[64, 64, 32],
        blocking_factor=[32, 32, 16],
        prob_relative=[3.0],
        fft=False,
    ),
    dict(
        id="Gauss2p5D",
        space_charge="Gauss2p5D",
        n_cell=[256, 256, 1],
        blocking_factor=[128, 128, 1],
        prob_relative=[3.0],
        fft=False,
    ),
]


@pytest.fixture(scope="function")
def sim(request):
    case = request.param

    # FFT-based solvers (IGF / 2D / 2p5D) only exist in an FFT-enabled build.
    # Probe this statically via Config so we can skip without a warm-up track.
    if case["fft"] and not Config.have_fft:
        pytest.skip("ImpactX built without FFT support; skipping FFT-only solver")

    class SimContextManager:
        def __enter__(self):
            self.sim = ImpactX()

            # set numerical parameters and IO control
            self.sim.particle_shape = 2  # B-spline order
            self.sim.diagnostics = False  # benchmarking
            self.sim.slice_step_diagnostics = False
            self.sim.tiny_profiler = False

            # space-charge solver configuration (must precede init_grids)
            self.sim.n_cell = case["n_cell"]
            self.sim.space_charge = case["space_charge"]
            if "poisson_solver" in case:
                self.sim.poisson_solver = case["poisson_solver"]
            self.sim.prob_relative = case["prob_relative"]
            if "blocking_factor" in case:
                bf_x, bf_y, bf_z = case["blocking_factor"]
                self.sim.blocking_factor_x = [bf_x]
                self.sim.blocking_factor_y = [bf_y]
                self.sim.blocking_factor_z = [bf_z]

            self.sim.init_grids()

            #   reference particle
            ref = self.sim.beam.ref
            ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

            #   particle bunch (Waterbag, matching the space-charge tests)
            distr = distribution.Waterbag(
                lambdaX=1.2154443728379865788e-3,
                lambdaY=1.2154443728379865788e-3,
                lambdaT=4.0956844276541331005e-4,
                lambdaPx=8.2274435782286157175e-4,
                lambdaPy=8.2274435782286157175e-4,
                lambdaPt=2.4415943602685364584e-3,
            )
            self.sim.add_particles(bunch_charge_C, distr, npart)
            assert self.sim.beam.total_number_of_particles() == npart

            # a single short CF cell, three slices -> three space-charge solves
            self.sim.lattice.extend(
                [
                    elements.ConstF(
                        name="cf1", ds=2.0, kx=1.0, ky=1.0, kt=1.0, nslice=3
                    ),
                ]
            )

            self.sim.backup_beam = self.sim.beam.make_alike()
            self.sim.backup_beam.arena = self.sim.beam.arena
            assert self.sim.backup_beam.total_number_of_particles() == 0

            self.sim.backup_beam.add_particles(self.sim.beam, local=True)
            assert self.sim.backup_beam.total_number_of_particles() == npart

            return self.sim

        def __exit__(self, exc_type, exc_value, traceback):
            # Work-around for https://github.com/AMReX-Codes/amrex/pull/5270
            self.sim.backup_beam.clear_particles()

            self.sim.finalize()
            del self.sim

    with SimContextManager() as sim:
        yield sim


def sc_setup(sim):
    """Fresh beam for each call of benchmark."""

    assert sim.backup_beam.total_number_of_particles() == npart

    # instead of drawing from the distribution again, restore the same initial
    # particles from the backup container.
    beam = sim.beam
    beam.clear_particles()
    beam.add_particles(sim.backup_beam, local=True)

    assert beam.total_number_of_particles() == npart

    # track_particles() takes no arguments
    return (), {}


@pytest.mark.parametrize(
    "sim", SC_CASES, indirect=True, ids=[c["id"] for c in SC_CASES]
)
def test_space_charge_solver(benchmark, sim):
    benchmark.pedantic(sim.track_particles, setup=partial(sc_setup, sim), rounds=rounds)
