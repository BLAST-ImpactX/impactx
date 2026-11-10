#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from math import factorial

import numpy as np
import pytest

from impactx import Config, ImpactX, elements


@pytest.mark.parametrize("skew", [False, True])
@pytest.mark.parametrize("m", [1, 2, 3, 4])
def test_exact_cfbend_multipole_symplectic(m, skew):
    """The map of a bend with a curvilinear multipole is symplectic: M^T J M = J.

    The Jacobian M is taken by central differences around a particle far off axis,
    where the multipole of index m (quadrupole: m = 1) has the same field at r = 5 cm
    for all orders.
    """
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    beam = sim.beam
    beam.ref.set_species("proton").set_kin_energy_MeV(2.0e3)

    is_single = Config.precision == "SINGLE"
    d = 1.0e-3 if is_single else 1.0e-6
    center = np.array([0.05, 0.01, -0.03, 0.005, 0.0, 0.0])  # x, px, y, py, t, pt
    phase_space = np.vstack([center + d * np.eye(6), center - d * np.eye(6)])
    x, px, y, py, t, pt = phase_space.T
    beam.add_n_particles(x, y, t, px, py, pt, beam.ref.qm_ratio_SI, 1.0e-12)

    k = [0.0] * 5
    k[m] = factorial(m) * 5.0e-3 / 0.05**m
    zero = [0.0] * 5
    elements.ExactCFbend(
        ds=0.15,
        k_normal=[2.0] + (zero if skew else k)[1:],
        k_skew=k if skew else zero,
    ).push(beam)

    df = beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
    out = df[
        [
            "position_x",
            "momentum_x",
            "position_y",
            "momentum_y",
            "position_t",
            "momentum_t",
        ]
    ].to_numpy(dtype=np.float64)

    sim.finalize()

    M = (out[:6] - out[6:]).T / (2.0 * d)
    J = np.kron(np.eye(3), np.array([[0.0, 1.0], [-1.0, 0.0]]))

    atol = 1.0e-8
    if is_single:
        # FIXME: the bend map and the multipole kicks cancel terms of order one in
        # single precision, M^T J M - J is up to 4e-4
        # https://github.com/BLAST-ImpactX/impactx/issues/1771
        atol = 2.0e-3
    np.testing.assert_allclose(M.T @ J @ M, J, rtol=0.0, atol=atol)
