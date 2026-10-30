#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# This set of tests is used for performance benchmarking of the beam
# distribution samplers (as micro-benchmarks). We use this file to rapidly
# evaluate performance changes when tuning the phase space and spin
# distributions in src/particles/distribution/ on CPUs and GPUs.
#
# Each round benchmarks ``sim.add_particles()`` into an emptied particle
# container: the distribution's ``initialize()``, the per-particle sampling
# kernel(s), and the copy of the sampled arrays into the particle container.
# Space charge is disabled, so no mesh resize or particle redistribution is
# measured.
#
# -*- coding: utf-8 -*-

import os
from functools import partial

import pytest

from impactx import ImpactX, distribution, twiss

# benchmark config
if os.environ.get("IS_CODESPEED_CPU_SIMULATION") == "1":
    # https://codspeed.io/docs/instruments/cpu/index
    rounds = 1
    npart = 10_000
else:
    rounds = 5
    npart = 1_000_000  # increase this to >10M or even 100M to avoid L1/L2/L3 cache effects on some hardware.

# a 1 GeV electron beam with an initial unnormalized rms emittance of 1 nm
kin_energy_MeV = 1.0e3  # reference energy
bunch_charge_C = 1.0e-9

# Phase space ellipse of all Twiss-parametrized distributions. The non-zero
# alphas keep the position-momentum correlation terms non-trivial.
beam_twiss = dict(
    beta_x=2.0,
    beta_y=1.5,
    beta_t=1.0,
    emitt_x=1.0e-09,
    emitt_y=1.0e-09,
    emitt_t=1.0e-06,
    alpha_x=-0.5,
    alpha_y=0.5,
    alpha_t=0.0,
)

# Thermal and bithermal beams of the EPAC2004 benchmarks
# (examples/epac2004_benchmarks/), which define their own reference particle
# and bunch charge.
thermal_kin_energy_MeV = 0.1
thermal_bunch_charge_C = 1.4285714285714285714e-10
THERMAL_CASES = {
    "thermal": dict(
        k=6.283185307179586,
        kT=36.0e-6,
        kT_halo=36.0e-6,
        normalize=0.41604661,
        normalize_halo=0.0,
        halo=0.0,
    ),
    "bithermal": dict(
        k=6.283185307179586,
        kT=36.0e-6,
        kT_halo=900.0e-6,
        normalize=0.54226,
        normalize_halo=0.08195,
        halo=0.05,
    ),
}


@pytest.fixture(scope="function")
def sim(request):
    """ImpactX sim with a 1 GeV electron reference particle and no particles.

    Parametrized indirectly by ``spin`` (default: off).
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

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

    yield sim

    sim.finalize()


def add_setup(sim, bunch_charge, distr, spin_distr):
    """Empty particle container for each call of benchmark."""

    beam = sim.beam
    beam.clear_particles()
    assert beam.total_number_of_particles() == 0

    return (bunch_charge, distr, npart, spin_distr), {}


def run_benchmark(benchmark, sim, distr, spin_distr=None, bunch_charge=bunch_charge_C):
    benchmark.pedantic(
        sim.add_particles,
        setup=partial(add_setup, sim, bunch_charge, distr, spin_distr),
        rounds=rounds,
    )
    assert sim.beam.total_number_of_particles() == npart


def test_Empty(benchmark, sim):
    """Baseline without sampling: allocation and copy into the particle container."""
    run_benchmark(benchmark, sim, distribution.Empty())


@pytest.mark.parametrize("cut", [0.0, 3.0], ids=["nocut", "cut"])
def test_Gaussian(benchmark, sim, cut):
    """3 Box-Muller pairs; ``cut`` samples a truncated radius per plane."""
    distr = distribution.Gaussian(**twiss(**beam_twiss), cutX=cut, cutY=cut, cutT=cut)
    run_benchmark(benchmark, sim, distr)


def test_Kurth4D(benchmark, sim):
    distr = distribution.Kurth4D(**twiss(**beam_twiss))
    run_benchmark(benchmark, sim, distr)


def test_Kurth6D(benchmark, sim):
    distr = distribution.Kurth6D(**twiss(**beam_twiss))
    run_benchmark(benchmark, sim, distr)


def test_KVdist(benchmark, sim):
    distr = distribution.KVdist(**twiss(**beam_twiss))
    run_benchmark(benchmark, sim, distr)


def test_Semigaussian(benchmark, sim):
    distr = distribution.Semigaussian(**twiss(**beam_twiss))
    run_benchmark(benchmark, sim, distr)


@pytest.mark.parametrize("params", THERMAL_CASES.values(), ids=THERMAL_CASES.keys())
def test_Thermal(benchmark, sim, params):
    """Includes the host-side radial CDF integration in ``initialize()`` and a
    binary search of that CDF per particle; ``bithermal`` mixes in a halo."""
    sim.beam.ref.set_species("proton").set_kin_energy_MeV(thermal_kin_energy_MeV)
    distr = distribution.Thermal(**params)
    run_benchmark(benchmark, sim, distr, bunch_charge=thermal_bunch_charge_C)


def test_Triangle(benchmark, sim):
    distr = distribution.Triangle(**twiss(**beam_twiss))
    run_benchmark(benchmark, sim, distr)


def test_Waterbag(benchmark, sim):
    distr = distribution.Waterbag(**twiss(**beam_twiss))
    run_benchmark(benchmark, sim, distr)


@pytest.mark.parametrize("sim", [True], indirect=True, ids=["spin"])
def test_SpinvMF(benchmark, sim):
    """Waterbag phase space plus von Mises-Fisher spin sampling: the spin
    sampling cost is the difference to ``test_Waterbag``."""
    distr = distribution.Waterbag(**twiss(**beam_twiss))
    spin_distr = distribution.SpinvMF(0.4, 0.9, 0.1)
    run_benchmark(benchmark, sim, distr, spin_distr=spin_distr)
