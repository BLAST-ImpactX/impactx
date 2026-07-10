#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# This set of tests is used for performance benchmarking of the reduced beam
# diagnostics as CodSpeed micro-benchmarks. We use these to rapidly evaluate
# performance changes when tuning the diagnostics reductions (see
# https://github.com/BLAST-ImpactX/impactx/issues/1102).
#
# The diagnostics calls benchmarked here are read-only with respect to the
# particle data, so no per-round particle reset is required. To benchmark an
# additional diagnostics call, add another ``test_*`` function that reuses the
# ``sim`` fixture below.
#
# -*- coding: utf-8 -*-

import os

import pytest

from impactx import ImpactX, distribution, twiss

# benchmark config
if os.environ.get("IS_CODESPEED_CPU_SIMULATION") == "1":
    # https://codspeed.io/docs/instruments/cpu/index
    rounds = 1
    npart = 10_000
else:
    rounds = 5
    npart = 1_000_000


@pytest.fixture(scope="function")
def sim(request):
    """ImpactX sim holding a 1 GeV electron Waterbag beam.

    Parametrized indirectly by ``spin`` (default: off). Diagnostics and space
    charge are disabled so that only the explicit diagnostics call under test
    is measured.
    """
    spin = getattr(request, "param", False)  # default to False if not parametrized

    sim = ImpactX()

    # numerical parameters and IO control
    sim.particle_shape = 2  # B-spline order
    sim.space_charge = False
    sim.diagnostics = False  # benchmarking
    sim.slice_step_diagnostics = False
    sim.spin = spin

    sim.init_grids()

    # a 1 GeV electron beam with an initial unnormalized rms emittance of 1 nm
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    distr = distribution.Waterbag(
        **twiss(
            beta_x=1.0,
            beta_y=1.0,
            beta_t=1.0,
            emitt_x=1.0e-09,
            emitt_y=1.0e-09,
            emitt_t=1.0e-06,
            alpha_x=0.0,
            alpha_y=0.0,
            alpha_t=0.0,
        )
    )
    bunch_charge_C = 1.0e-9
    if spin:
        spin_vectors = distribution.SpinvMF(0.4, 0.9, 0.1)
        sim.add_particles(bunch_charge_C, distr, npart, spin_vectors)
    else:
        sim.add_particles(bunch_charge_C, distr, npart)
    assert sim.beam.total_number_of_particles() == npart

    yield sim

    sim.finalize()


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
@pytest.mark.parametrize("eigenemittances", [False, True], ids=["noeig", "eig"])
def test_beam_moments(benchmark, sim, eigenemittances):
    """Benchmark the reduced beam characteristics (means, sigmas, min/max,
    emittances, Twiss parameters, dispersion) from the particle distribution.

    - spin / nospin: the three spin components are currently reduced regardless
      of whether spin tracking is enabled.
    - noeig / eig: ``eig`` additionally computes the three eigenemittances from
      the full 6x6 beam covariance matrix (a per-call dense eigen-solve).
    """
    sim.eigenemittances = eigenemittances
    benchmark.pedantic(sim.beam.beam_moments, rounds=rounds)
