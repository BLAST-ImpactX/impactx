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

single_precision = Config.precision_particles == "SINGLE"

# One rounding error in a radicand, or one draw at the end of its range, is
# enough for a non-finite particle. Rounding errors in a radicand are on the
# order of 1e-7 per particle in single precision.
npart_rounding = 20_000_000 if single_precision else 1_000_000
# A draw at the end of its range has P = 2.8e-8 in single precision on GPU. On
# CPU, it has P = 2^-32, which is out of reach of this test.
npart_draw = 20_000_000 if single_precision and Config.have_gpu else 1_000_000

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
    # ... and the fraction outside the cut underflows to 0
    "Gaussian_cut50": lambda: gaussian_cut(50.0),
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


def radius(df, plane, sigma, sigma_p):
    """Radius of each particle in a phase plane, in units of the sampled sigma."""
    position = df[f"position_{plane}"].to_numpy(dtype=np.float64) / sigma
    momentum = df[f"momentum_{plane}"].to_numpy(dtype=np.float64) / sigma_p
    return np.sqrt(position**2 + momentum**2)


@pytest.mark.parametrize("name", DISTRIBUTIONS)
def test_distribution_finite(sim, name):
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    npart = npart_rounding if name in ("Kurth4D", "Kurth6D") else npart_draw
    sim.add_particles(1.0e-9, DISTRIBUTIONS[name](), npart)

    assert_finite_moments(sim)


@pytest.mark.parametrize("cut", [1.0, 3.0])
def test_gaussian_cut_radius(sim, cut):
    """The radius in a phase plane ends at the cut, in units of the untruncated sigma."""
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    sigma = {"x": 2.0e-3, "px": 3.0e-3}
    distr = distribution.Gaussian(
        lambdaX=sigma["x"],
        lambdaY=1.0e-3,
        lambdaT=1.0e-3,
        lambdaPx=sigma["px"],
        lambdaPy=1.0e-3,
        lambdaPt=1.0e-3,
        cutX=cut,
    )
    sim.add_particles(1.0e-9, distr, 1_000_000)

    df = sim.beam.to_df(local=True)

    # the sampled plane is rescaled to unit variance, by the variance of the truncation
    f = cut**2 / 2.0
    r_cut = cut / np.sqrt(1.0 - f / np.expm1(f))

    r_max = radius(df, "x", sigma["x"], sigma["px"]).max()
    assert r_max <= r_cut * (1.0 + 1.0e-5)
    assert r_max >= r_cut * (1.0 - 1.0e-2)


@pytest.mark.skipif(
    not single_precision, reason="the rounding is out of reach in double precision"
)
def test_gaussian_cut_tail(sim):
    """The tail of a large cut is resolved below the rounding of 1.

    The radius is r^2 = -2 ln(w). A w computed as 1 - u*x is a multiple of 2^-24
    in the tail, in single precision, which gives discrete radii that end at
    5.77 sigma.
    """
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    sigma = 1.0e-3
    distr = distribution.Gaussian(
        lambdaX=sigma,
        lambdaY=sigma,
        lambdaT=sigma,
        lambdaPx=sigma,
        lambdaPy=sigma,
        lambdaPt=sigma,
        cutX=7.0,
        cutY=7.0,
        cutT=7.0,
    )
    sim.add_particles(1.0e-9, distr, 4_000_000)

    df = sim.beam.to_df(local=True)
    r = np.concatenate([radius(df, plane, sigma, sigma) for plane in "xyt"])

    # w in units of 2^-24; this cut does not rescale the variance in single precision
    k = np.exp(-(r**2) / 2.0) / np.finfo(np.float32).epsneg
    # on average 46 particles
    tail = k[k < 64.0]
    assert tail.size > 10
    # about 2% of a continuous w are this close to a multiple of 2^-24
    assert np.mean(np.abs(tail - np.rint(tail)) < 0.01) < 0.5


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
    sim.add_particles(1.4285714285714285714e-10, thermal, npart_draw)

    assert_finite_moments(sim)
