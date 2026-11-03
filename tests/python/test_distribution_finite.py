#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, distribution, twiss

# One draw at the end of its range, or one rounding error in a radicand, is
# enough for a non-finite particle. These are rare per particle, in single
# precision on the order of 1e-7.
npart = 20_000_000 if Config.precision == "SINGLE" else 1_000_000

beam_twiss = twiss(
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


def gaussian_cut(cut):
    return distribution.Gaussian(**beam_twiss, cutX=cut, cutY=cut, cutT=cut)


DISTRIBUTIONS = {
    "Gaussian": lambda: distribution.Gaussian(**beam_twiss),
    "Gaussian_cut3": lambda: gaussian_cut(3.0),
    # the fraction of the beam inside the cut rounds to exactly 1 in single precision
    "Gaussian_cut6": lambda: gaussian_cut(6.0),
    # ... and in double precision
    "Gaussian_cut10": lambda: gaussian_cut(10.0),
    "Kurth4D": lambda: distribution.Kurth4D(**beam_twiss),
    "Kurth6D": lambda: distribution.Kurth6D(**beam_twiss),
    "KVdist": lambda: distribution.KVdist(**beam_twiss),
    "Semigaussian": lambda: distribution.Semigaussian(**beam_twiss),
    "Triangle": lambda: distribution.Triangle(**beam_twiss),
    "Waterbag": lambda: distribution.Waterbag(**beam_twiss),
}


@pytest.fixture
def sim():
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    yield sim

    sim.finalize()


def assert_finite_moments(sim):
    """One non-finite particle makes the sum over all particles non-finite.

    The means are plain sums. The RMS sizes are not suitable here: they are
    clamped with a comparison, which turns a NaN to 0 with fast-math.
    """
    moments = sim.beam.beam_moments()
    for moment in ["mean_x", "mean_y", "mean_t", "mean_px", "mean_py", "mean_pt"]:
        assert np.isfinite(moments[moment]), f"{moment} = {moments[moment]}"


@pytest.mark.parametrize("name", DISTRIBUTIONS)
def test_distribution_finite(sim, name):
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    sim.add_particles(1.0e-9, DISTRIBUTIONS[name](), npart)

    assert_finite_moments(sim)


def test_thermal_finite(sim):
    sim.beam.ref.set_species("proton").set_kin_energy_MeV(0.1)
    thermal = distribution.Thermal(
        k=6.283185307179586,
        kT=36.0e-6,
        kT_halo=36.0e-6,
        normalize=0.41604661,
        normalize_halo=0.0,
        halo=0.0,
    )
    sim.add_particles(1.4285714285714285714e-10, thermal, npart)

    assert_finite_moments(sim)
