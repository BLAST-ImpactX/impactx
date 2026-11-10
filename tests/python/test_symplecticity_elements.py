#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# This set of tests verifies that the particle map of each element is
# symplectic: its Jacobian M satisfies M^T J M = J.
#
# -*- coding: utf-8 -*-

from math import factorial

import numpy as np
import pytest

from impactx import Config, ImpactX, elements

# canonical pairs (x, px), (y, py), (t, pt)
PHASE_COLS = [
    "position_x",
    "momentum_x",
    "position_y",
    "momentum_y",
    "position_t",
    "momentum_t",
]
J = np.kron(np.eye(3), np.array([[0.0, 1.0], [-1.0, 0.0]]))

# a particle far off axis, where nonlinear fields are strong: x, px, y, py, t, pt
CENTER = np.array([0.05, 0.01, -0.03, 0.005, 0.0, 0.0])

# finite-difference step of the Jacobian, and tolerance of M^T J M - J
if Config.precision == "SINGLE":
    step = 1.0e-3
    atol = 3.0e-4
else:
    step = 1.0e-4
    atol = 2.0e-9


@pytest.fixture(scope="function")
def sim():
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("proton").set_kin_energy_MeV(2.0e3)

    yield sim

    sim.finalize()


def jacobian(el, sim, center=CENTER):
    """Jacobian of the particle map of an element, by fourth-order central differences"""
    offsets = [2.0, 1.0, -1.0, -2.0]
    phase_space = np.vstack([center + o * step * np.eye(6) for o in offsets])
    x, px, y, py, t, pt = phase_space.T

    beam = sim.beam
    beam.add_n_particles(x, y, t, px, py, pt, beam.ref.qm_ratio_SI, 1.0e-12)
    el.push(beam)

    df = beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
    plus2, plus1, minus1, minus2 = np.split(
        df[PHASE_COLS].to_numpy(dtype=np.float64), len(offsets)
    )
    return (-plus2 + 8.0 * plus1 - 8.0 * minus1 + minus2).T / (12.0 * step)


def check_symplectic(el, sim, center=CENTER, atol=atol):
    """Assert that the particle map of an element is symplectic at a point"""
    M = jacobian(el, sim, center)
    np.testing.assert_allclose(M.T @ J @ M, J, rtol=0.0, atol=atol)


@pytest.mark.parametrize("skew", [False, True], ids=["normal", "skew"])
@pytest.mark.parametrize("m", range(1, 12))
def test_ExactCFbend(sim, m, skew):
    """Each curvilinear multipole of index m (quadrupole: m = 1), with the same
    field at r = 5 cm for all orders
    """
    k = [0.0] * (m + 1)
    k[m] = factorial(m) * 5.0e-3 / 0.05**m
    zero = [0.0] * (m + 1)
    el = elements.ExactCFbend(
        ds=0.15,
        k_normal=[2.0] + (zero if skew else k)[1:],
        k_skew=k if skew else zero,
    )

    tol = atol
    if Config.precision == "SINGLE":
        # FIXME: the bend map cancels terms of order rc in single precision,
        # M^T J M - J is up to 5e-4
        # https://github.com/BLAST-ImpactX/impactx/issues/1771
        tol = 2.0e-3
    check_symplectic(el, sim, atol=tol)
