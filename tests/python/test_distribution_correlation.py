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

from impactx import ImpactX, create_envelope, distribution, twiss

SAMPLERS = {
    "Gaussian": distribution.Gaussian,
    "Kurth4D": distribution.Kurth4D,
    "Kurth6D": distribution.Kurth6D,
    "KVdist": distribution.KVdist,
    "Semigaussian": distribution.Semigaussian,
    "Triangle": distribution.Triangle,
    "Waterbag": distribution.Waterbag,
}

# 1 - mu^2 is 7.7e-8, or 2e-8, while a float32 rounds mu to a distance from 1 of 6e-8 or 0
CORRELATIONS = [0.999999961419755, 0.99999999]

NPART = 200_000
# RMS sizes and emittances from random sampling, for this number of particles
RTOL = 5.0 / np.sqrt(NPART)

LAMBDA_T = 5.557e-9
LAMBDA_PT = 7.2e-2


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


def longitudinal_moments(df):
    """RMS sizes and emittance in (t,pt), in double precision from the particles."""
    t = df["position_t"].to_numpy(dtype=np.float64)
    pt = df["momentum_t"].to_numpy(dtype=np.float64)
    t = t - t.mean()
    pt = pt - pt.mean()
    slope = np.dot(t, pt) / np.dot(t, t)
    sigma_t = np.sqrt(np.mean(t * t))
    sigma_pt = np.sqrt(np.mean(pt * pt))
    emittance_t = sigma_t * np.sqrt(np.mean((pt - slope * t) ** 2))
    return sigma_t, sigma_pt, emittance_t


def design_moments(mu):
    root = np.sqrt((1.0 - mu) * (1.0 + mu))
    return LAMBDA_T / root, LAMBDA_PT / root, LAMBDA_T * LAMBDA_PT / root


def beam_ellipse(mu):
    return dict(
        lambdaX=1.0e-5,
        lambdaY=1.0e-5,
        lambdaT=LAMBDA_T,
        lambdaPx=1.0e-5,
        lambdaPy=1.0e-5,
        lambdaPt=LAMBDA_PT,
        mutpt=mu,
    )


@pytest.mark.parametrize("mu", CORRELATIONS)
@pytest.mark.parametrize("name", SAMPLERS)
def test_near_degenerate_beam(sim, name, mu):
    """The sampled beam has the requested sizes, also for a correlation close to 1."""
    sim.add_particles(1.0e-9, SAMPLERS[name](**beam_ellipse(mu)), NPART)

    moments = longitudinal_moments(sim.beam.to_df(local=True))

    np.testing.assert_allclose(moments, design_moments(mu), rtol=RTOL)


def test_near_degenerate_beam_negative(sim):
    mu = 0.999999961419755
    sim.add_particles(1.0e-9, distribution.Gaussian(**beam_ellipse(-mu)), NPART)

    moments = longitudinal_moments(sim.beam.to_df(local=True))

    np.testing.assert_allclose(np.abs(moments), design_moments(mu), rtol=RTOL)


@pytest.mark.parametrize("alpha", [3.6e3, 1.0e4, 1.0e5])
def test_near_degenerate_beam_twiss(sim, alpha):
    """A large alpha in the Courant-Snyder inputs gives a correlation close to 1."""
    emittance, beta = 1.0e-9, 1.0
    beam_twiss = twiss(
        beta_x=2.0,
        beta_y=2.0,
        beta_t=beta,
        emitt_x=1.0e-9,
        emitt_y=1.0e-9,
        emitt_t=emittance,
        alpha_x=0.0,
        alpha_y=0.0,
        alpha_t=alpha,
    )
    sim.add_particles(1.0e-9, distribution.Gaussian(**beam_twiss), NPART)

    moments = longitudinal_moments(sim.beam.to_df(local=True))

    gamma = (1.0 + alpha**2) / beta
    design = (np.sqrt(emittance * beta), np.sqrt(emittance * gamma), emittance)
    np.testing.assert_allclose(moments, design, rtol=RTOL)


@pytest.mark.parametrize("mu", CORRELATIONS)
def test_near_degenerate_envelope(sim, mu):
    """The covariance matrix of an envelope has the requested sizes."""
    envelope = create_envelope(distribution.Gaussian(**beam_ellipse(mu)), 1.0e-9)
    moments = envelope.beam_moments(sim.beam.ref)

    sigma_t, sigma_pt, _ = design_moments(mu)
    np.testing.assert_allclose(
        [moments["sigma_t"], moments["sigma_pt"]], [sigma_t, sigma_pt], rtol=1.0e-5
    )


@pytest.mark.parametrize("mu", [1.0, -1.0, 1.5])
@pytest.mark.parametrize("plane", ["muxpx", "muypy", "mutpt"])
def test_invalid_correlation(plane, mu):
    parameters = beam_ellipse(0.0)
    parameters[plane] = mu
    with pytest.raises(RuntimeError, match="correlation"):
        distribution.Gaussian(**parameters)
