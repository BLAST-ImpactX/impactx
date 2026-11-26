#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""The covariance matrix of a particle beam."""

import math

import numpy as np
import pytest

from impactx import Config, ImpactX, distribution, elements

COLUMNS = (
    "position_x",
    "momentum_x",
    "position_y",
    "momentum_y",
    "position_t",
    "momentum_t",
)


def make_sim(npart):
    """A 100 MeV electron beam with correlations in each plane."""
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(100.0)
    distr = distribution.Waterbag(
        lambdaX=1.0e-4,
        lambdaY=2.0e-4,
        lambdaT=1.0e-3,
        lambdaPx=3.0e-5,
        lambdaPy=1.0e-5,
        lambdaPt=2.0e-3,
        muxpx=-0.5,
        muypy=0.3,
        mutpt=0.4,
    )
    if npart > 0:
        sim.add_particles(1.0e-9, distr, npart)
    return sim


def to_numpy(cm):
    """A Map6x6 (indices start at 1) as a numpy array."""
    return np.array([[cm[i, j] for j in range(1, 7)] for i in range(1, 7)])


def weighted_covariance(df):
    """Weighted covariance of x, px, y, py, t, pt, in two passes."""
    w = df["weighting"].to_numpy(dtype=np.float64)
    data = np.vstack([df[c].to_numpy(dtype=np.float64) for c in COLUMNS])
    mean = (data * w).sum(axis=1) / w.sum()
    d = data - mean[:, None]
    return (d * w) @ d.T / w.sum()


def test_covariance_matrix_coupled_beam():
    """All 21 moments, including the couplings between the planes"""
    sim = make_sim(10000)
    # couple the planes: a solenoid couples x and y, a sector bend x and t
    elements.Sol(ds=0.5, ks=2.0).push(sim.beam)
    elements.Sbend(ds=0.5, rc=2.0, nslice=4).push(sim.beam)

    cm = to_numpy(sim.beam.covariance_matrix())
    np.testing.assert_array_equal(cm, cm.T)

    # the same moments as the reduced beam characteristics
    rbc = sim.beam.beam_moments()
    rtol = 1.0e-5 if Config.precision == "SINGLE" else 1.0e-12
    for i, name in enumerate(("x", "px", "y", "py", "t", "pt")):
        assert math.isclose(cm[i, i], rbc[f"sigma_{name}"] ** 2, rel_tol=rtol)

    # all particles, gathered on MPI rank 0
    df = sim.beam.to_df(local=False)
    if df is not None:
        exact = weighted_covariance(df)
        scale = np.sqrt(np.outer(np.diag(exact), np.diag(exact)))
        # the planes are coupled, beyond the dispersion
        for i, j in ((0, 2), (1, 3), (0, 4), (1, 4)):
            assert abs(exact[i, j]) > 0.01 * scale[i, j]
        # each moment relative to its scale sqrt(<u^2><v^2>)
        atol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-10
        np.testing.assert_allclose(cm / scale, exact / scale, rtol=0.0, atol=atol)

    sim.finalize()


def test_covariance_matrix_without_beam():
    """Without particles, the covariance matrix is undefined"""
    sim = make_sim(0)
    with pytest.raises(RuntimeError, match="no particle weight"):
        sim.beam.covariance_matrix()
    sim.finalize()


def test_covariance_matrix_without_reference():
    """Without a reference particle energy, the covariance matrix is undefined"""
    sim = make_sim(100)
    sim.beam.ref.pt = 0.0  # gamma = 0: not initialized
    with pytest.raises(RuntimeError, match="reference particle is not initialized"):
        sim.beam.covariance_matrix()
    sim.finalize()
