#!/usr/bin/env python3
#
# Copyright 2026 The ImpactX Community
#
# Authors: OpenAI Codex
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, RefPart, create_envelope, distribution, elements


def phase_tolerances():
    if Config.precision == "SINGLE":
        return 1.0e-7, 5.0e-5
    return 1.0e-14, 1.0e-8


def test_create_envelope_initializes_centroid_from_distribution_means():
    atol, rtol = phase_tolerances()

    base_kwargs = dict(
        lambdaX=4.0e-5,
        lambdaY=5.0e-5,
        lambdaT=1.0e-3,
        lambdaPx=1.0e-5,
        lambdaPy=3.0e-5,
        lambdaPt=2.0e-3,
        muxpx=-0.2,
        muypy=0.3,
        mutpt=-0.1,
    )
    mean_kwargs = dict(
        meanX=1.5e-3,
        meanPx=-2.5e-4,
        meanY=-3.5e-3,
        meanPy=4.5e-4,
        meanT=5.5e-3,
        meanPt=-6.5e-4,
    )

    env = create_envelope(distribution.Gaussian(**base_kwargs, **mean_kwargs), None)
    env_zero_mean = create_envelope(distribution.Gaussian(**base_kwargs), None)

    assert np.allclose(
        np.asarray(env.centroid),
        np.array(
            [
                mean_kwargs["meanX"],
                mean_kwargs["meanPx"],
                mean_kwargs["meanY"],
                mean_kwargs["meanPy"],
                mean_kwargs["meanT"],
                mean_kwargs["meanPt"],
            ]
        ),
        atol=atol,
        rtol=rtol,
    )
    assert np.allclose(
        env.envelope.to_numpy(),
        env_zero_mean.envelope.to_numpy(),
        atol=atol,
        rtol=rtol,
    )


def test_create_envelope_still_rejects_nonzero_dispersion():
    with pytest.raises(RuntimeError, match="non-zero dispersion"):
        create_envelope(
            distribution.Gaussian(
                lambdaX=4.0e-5,
                lambdaY=5.0e-5,
                lambdaT=1.0e-3,
                lambdaPx=1.0e-5,
                lambdaPy=3.0e-5,
                lambdaPt=2.0e-3,
                dispX=1.0e-3,
            ),
            None,
        )


def test_track_envelope_centroid_propagates_constant_kicks():
    atol, rtol = phase_tolerances()

    sim = ImpactX()
    sim.init_grids()

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    sim.init_envelope(
        ref,
        distribution.Gaussian(
            lambdaX=4.0e-5,
            lambdaY=5.0e-5,
            lambdaT=1.0e-3,
            lambdaPx=1.0e-5,
            lambdaPy=3.0e-5,
            lambdaPt=2.0e-3,
        ),
    )

    k1 = 9.0
    dx = 2.0e-3
    drift_ds = 0.4
    k2 = 7.0

    sim.lattice.extend(
        [
            elements.Multipole(multipole=3, K_normal=k1, K_skew=0.0, dx=dx),
            elements.Drift(ds=drift_ds),
            elements.Multipole(multipole=3, K_normal=k2, K_skew=0.0),
        ]
    )

    sim.track_envelope()

    dpx_1 = -0.5 * k1 * dx * dx
    x_after_drift = drift_ds * dpx_1
    dpx_2 = -0.5 * k2 * x_after_drift * x_after_drift

    assert np.allclose(
        np.asarray(sim.envelope.centroid),
        np.array([x_after_drift, dpx_1 + dpx_2, 0.0, 0.0, 0.0, 0.0]),
        atol=atol,
        rtol=rtol,
    )

    sim.finalize()
