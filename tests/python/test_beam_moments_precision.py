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

from impactx import Config, ImpactX, distribution

single_precision = Config.precision_particles == "SINGLE"

# relative error of the moments in the sums over 1M particles, against the same moments in
# double precision: this is the rounding of the sums, not the sampling
RTOL = 1.0e-3 if single_precision else 1.0e-8

BUNCH_CHARGE_C = 1.0e-9
ELECTRON_QM_EEV = -1.0 / 0.510998950 / 1.0e6  # charge/mass in e / eV

# 1 - mu^2 = 7.7e-8 is the input of examples/chicane
NEAR_DEGENERATE = 0.999999961419755


@pytest.fixture
def sim():
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()
    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    yield sim

    sim.finalize()


def exact_moments(df, weight=None):
    """The beam moments in double precision from the particles, by a stable two-pass."""
    w = np.ones(len(df)) if weight is None else weight
    w = w / w.sum()

    coord = {k: df[f"position_{k}"].to_numpy(dtype=np.float64) for k in "xyt"}
    coord.update(
        {f"p{k}": df[f"momentum_{k}"].to_numpy(dtype=np.float64) for k in "xyt"}
    )

    m = {}
    for name, u in coord.items():
        m[f"mean_{name}"] = np.sum(w * u)
        coord[name] = u - m[f"mean_{name}"]
        m[f"sigma_{name}"] = np.sqrt(np.sum(w * coord[name] ** 2))

    def emittance(q, p):
        slope = np.sum(w * q * p) / np.sum(w * q * q)
        return np.sqrt(np.sum(w * q * q) * np.sum(w * (p - slope * q) ** 2))

    for k in "xyt":
        m[f"emittance_{k}"] = emittance(coord[k], coord[f"p{k}"])

    pt = coord["pt"]
    pt2 = np.sum(w * pt * pt)
    for name in ("x", "px", "y", "py"):
        m[f"dispersion_{name}"] = -np.sum(w * coord[name] * pt) / pt2

    # Courant-Snyder functions, with the dispersion removed in x and y
    for k in "xy":
        q = coord[k] + m[f"dispersion_{k}"] * pt
        p = coord[f"p{k}"] + m[f"dispersion_p{k}"] * pt
        eps = emittance(q, p)
        m[f"beta_{k}"] = np.sum(w * q * q) / eps
        m[f"alpha_{k}"] = -np.sum(w * q * p) / eps
    t = coord["t"]
    m["beta_t"] = np.sum(w * t * t) / m["emittance_t"]
    m["alpha_t"] = -np.sum(w * t * pt) / m["emittance_t"]

    return m, coord


def assert_moments_close(moments, exact, coord, rtol=RTOL):
    """Each moment, in the scale it is known to: a mean relative to the RMS size of its coordinate."""
    for key, value in exact.items():
        if key.startswith("mean_"):
            name = key.removeprefix("mean_")
            scale = np.sqrt(np.mean(coord[name] ** 2))
        elif key.startswith("dispersion_"):
            name = key.removeprefix("dispersion_")
            scale = np.sqrt(np.mean(coord[name] ** 2) / np.mean(coord["pt"] ** 2))
        elif key.startswith("alpha_"):
            # alpha of a plane without correlation is sampling noise around 0
            scale = max(1.0, abs(value))
        else:
            scale = abs(value)
        error = abs(moments[key] - value) / scale
        assert error <= rtol, f"{key} = {moments[key]}, exact {value}: {error}"


@pytest.mark.parametrize(
    "correlations",
    [
        dict(mutpt=NEAR_DEGENERATE),
        dict(muxpx=NEAR_DEGENERATE),
        dict(muypy=-NEAR_DEGENERATE),
        dict(muxpx=NEAR_DEGENERATE, muypy=-NEAR_DEGENERATE, mutpt=NEAR_DEGENERATE),
        dict(muxpx=0.9999, muypy=0.99, mutpt=0.5),
        dict(muxpx=0.5, muypy=-0.3),
    ],
)
def test_beam_moments_correlated(sim, correlations):
    """The moments of a phase plane with a strong correlation, which cancel in a plain sum."""
    distr = distribution.Gaussian(
        lambdaX=1.0e-5,
        lambdaY=2.0e-5,
        lambdaT=5.557e-9,
        lambdaPx=3.0e-5,
        lambdaPy=4.0e-5,
        lambdaPt=7.2e-5,
        meanX=1.0e-5,
        meanPx=-2.0e-5,
        meanT=3.0e-8,
        meanPt=1.0e-5,
        **correlations,
    )
    sim.add_particles(BUNCH_CHARGE_C, distr, 1_000_000)

    exact, coord = exact_moments(sim.beam.to_df(local=True))

    # double precision does not shear the momenta (see shear_momenta in
    # ReducedBeamCharacteristics.cpp): the moments of a near-degenerate plane keep the
    # rounding of the sums times 1/(1 - r^2), which is 2.4e-7 here
    near_degenerate = NEAR_DEGENERATE in map(abs, correlations.values())
    rtol = 1.0e-6 if near_degenerate and not single_precision else RTOL
    assert_moments_close(sim.beam.beam_moments(), exact, coord, rtol=rtol)


def test_beam_moments_dispersive(sim):
    """The moments of all planes and the correlations between them, recovered from a shear."""
    distr = distribution.Gaussian(
        lambdaX=1.0e-5,
        lambdaY=2.0e-5,
        lambdaT=1.0e-6,
        lambdaPx=3.0e-5,
        lambdaPy=4.0e-5,
        lambdaPt=1.0e-5,
        muxpx=0.999,
        muypy=-0.999,
        mutpt=0.99999,
        dispX=1.0e-3,
        dispPx=2.0e-4,
        dispY=-5.0e-4,
        dispPy=1.0e-4,
    )
    sim.add_particles(BUNCH_CHARGE_C, distr, 1_000_000)

    exact, coord = exact_moments(sim.beam.to_df(local=True))

    # the dispersion takes the largest part of the size in x and y, which cancels
    assert_moments_close(sim.beam.beam_moments(), exact, coord, rtol=30.0 * RTOL)


@pytest.mark.parametrize("npart", [1, 2, 3, 4, 31, 32, 33])
def test_beam_moments_few_particles(sim, npart):
    """The sample that the shifts are taken from can have any size."""
    rng = np.random.default_rng(1)
    x, y, t = 1.0e-3 * rng.normal(size=(3, npart))
    px, py, pt = 2.0e-3 * x + 1.0e-6 * rng.normal(size=(3, npart))
    sim.beam.add_n_particles(x, y, t, px, py, pt, ELECTRON_QM_EEV, BUNCH_CHARGE_C)

    moments = sim.beam.beam_moments()

    df = sim.beam.to_df(local=True)
    px_exact = df["momentum_x"].to_numpy(dtype=np.float64)
    assert moments["mean_px"] == pytest.approx(np.mean(px_exact), rel=1.0e-5)
    assert moments["sigma_px"] == pytest.approx(np.std(px_exact), rel=1.0e-2, abs=1e-9)


def test_beam_moments_thin_sample(sim):
    """The first particles in memory do not represent the beam, here all from a thin slice."""
    npart = 1_000_000
    rng = np.random.default_rng(1)
    t = 1.0e-3 * rng.normal(size=npart)
    # uncorrelated: a slope from a thin slice in t would add a correlation that is not there
    pt = 2.0e-3 * rng.normal(size=npart)
    x, y, px, py = 1.0e-3 * rng.normal(size=(4, npart))
    slice_first = np.argsort(np.abs(t))
    sim.beam.add_n_particles(
        x[slice_first],
        y[slice_first],
        t[slice_first],
        px[slice_first],
        py[slice_first],
        pt[slice_first],
        ELECTRON_QM_EEV,
        BUNCH_CHARGE_C,
    )

    exact, coord = exact_moments(sim.beam.to_df(local=True))

    assert_moments_close(sim.beam.beam_moments(), exact, coord)


def test_beam_moments_constant_coordinate(sim):
    """A plane without a spread has no size, also if it is not at the reference particle."""
    distr = distribution.Gaussian(
        lambdaX=1.0e-5,
        lambdaY=2.0e-5,
        lambdaT=1.0e-5,
        lambdaPx=3.0e-5,
        lambdaPy=4.0e-5,
        lambdaPt=0.0,
        meanPt=1.0e-3 / 3.0,
    )
    sim.add_particles(BUNCH_CHARGE_C, distr, 1000)

    moments = sim.beam.beam_moments()

    assert moments["sigma_pt"] == 0.0
    assert moments["emittance_t"] == 0.0
    assert moments["emittance_x"] > 0.0


@pytest.mark.parametrize("npart", [1_000_000, 4_000_000])
def test_beam_moments_charge(sim, npart):
    """Summing the same weight in a floating point number rounds in the same direction."""
    distr = distribution.Gaussian(
        lambdaX=1.0e-5,
        lambdaY=1.0e-5,
        lambdaT=1.0e-5,
        lambdaPx=1.0e-5,
        lambdaPy=1.0e-5,
        lambdaPt=1.0e-5,
    )
    sim.add_particles(BUNCH_CHARGE_C, distr, npart)

    moments = sim.beam.beam_moments()

    assert abs(moments["charge_C"]) == pytest.approx(
        BUNCH_CHARGE_C, rel=1.0e-6 if single_precision else 1.0e-13
    )


def test_beam_moments_weights(sim):
    """Particles of different weights."""
    npart = 200_000
    rng = np.random.default_rng(1)
    x, y, t, px, py, pt = 1.0e-3 * rng.normal(size=(6, npart))
    w = rng.uniform(0.5e4, 1.5e4, size=npart)
    sim.beam.add_n_particles(x, y, t, px, py, pt, ELECTRON_QM_EEV, w=w)

    df = sim.beam.to_df(local=True)
    exact, coord = exact_moments(df, df["weighting"].to_numpy(dtype=np.float64))

    moments = sim.beam.beam_moments()
    assert_moments_close(moments, exact, coord)
    elementary_charge_C = 1.602176634e-19
    assert abs(moments["charge_C"]) == pytest.approx(
        np.sum(w) * elementary_charge_C, rel=1.0e-4 if single_precision else 1.0e-13
    )


def test_beam_moments_heavy_first_particles(sim):
    """The first particles in memory, which the shifts are taken from, outweigh the others."""
    npart = 1_000_000
    rng = np.random.default_rng(1)
    x, y, t, px, py, pt = 1.0e-3 * rng.normal(size=(6, npart))
    w = np.ones(npart)
    w[:32] = 1.0e4
    sim.beam.add_n_particles(x, y, t, px, py, pt, ELECTRON_QM_EEV, w=w)

    df = sim.beam.to_df(local=True)
    exact, coord = exact_moments(df, df["weighting"].to_numpy(dtype=np.float64))

    moments = sim.beam.beam_moments()
    assert_moments_close(moments, exact, coord)
    elementary_charge_C = 1.602176634e-19
    assert abs(moments["charge_C"]) == pytest.approx(
        np.sum(w) * elementary_charge_C, rel=1.0e-6 if single_precision else 1.0e-13
    )


@pytest.mark.skipif(
    single_precision,
    # FIXME: the coefficients of the characteristic polynomial cancel in single precision
    reason="https://github.com/BLAST-ImpactX/impactx/issues/1481",
)
def test_eigenemittances_near_degenerate(sim):
    """Without a coupling, the eigenemittances are the emittances of the planes."""
    sim.eigenemittances = True
    distr = distribution.Gaussian(
        lambdaX=1.0e-5,
        lambdaY=2.0e-5,
        lambdaT=5.557e-9,
        lambdaPx=3.0e-5,
        lambdaPy=4.0e-5,
        lambdaPt=7.2e-5,
        muxpx=0.99,
        muypy=-0.99,
        mutpt=NEAR_DEGENERATE,
    )
    sim.add_particles(BUNCH_CHARGE_C, distr, 1_000_000)

    moments = sim.beam.beam_moments()

    # statistical correlations of 1/sqrt(N) between the planes couple them
    np.testing.assert_allclose(
        sorted(
            [moments["emittance_1"], moments["emittance_2"], moments["emittance_3"]]
        ),
        sorted(
            [moments["emittance_xn"], moments["emittance_yn"], moments["emittance_tn"]]
        ),
        rtol=5.0e-3,
    )
