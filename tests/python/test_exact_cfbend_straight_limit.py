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

PHASE_SPACE_COLUMNS = ["position_x", "momentum_x", "position_y", "momentum_y"]


def track(element, radius):
    """Phase space of a few particles within a radius after the element"""
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    beam = sim.beam
    beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    x = radius * np.array([0.75, -0.5, 1.0, -0.25])
    y = radius * np.array([0.25, -0.75, -0.5, 1.0])
    px = np.array([1.0e-4, -2.0e-4, 0.0, 3.0e-4])
    py = np.array([-1.0e-4, 0.0, 2.0e-4, 1.0e-4])
    zero = np.zeros_like(x)
    beam.add_n_particles(x, y, zero, px, py, zero, beam.ref.qm_ratio_SI, 1.0e-12)

    element.push(beam)

    df = beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
    phase_space = df[PHASE_SPACE_COLUMNS].to_numpy(dtype=np.float64)
    sim.finalize()
    return phase_space


def test_exact_cfbend_zero_curvature():
    """Without curvature, the multipoles of all orders kick as in ExactMultipole."""
    k_normal = [0.0] + [factorial(m) * 1.0e-3 / 5.0e-3**m for m in range(1, 12)]
    k_skew = [0.0] + [
        (-1) ** m * factorial(m) * 0.5e-3 / 5.0e-3**m for m in range(1, 12)
    ]
    kwargs = dict(ds=0.1, k_normal=k_normal, k_skew=k_skew, int_order=2, mapsteps=10)

    bend = track(elements.ExactCFbend(**kwargs), 4.0e-3)
    straight = track(elements.ExactMultipole(**kwargs), 4.0e-3)

    rtol = 1.0e-5 if Config.precision == "SINGLE" else 1.0e-12
    np.testing.assert_allclose(bend, straight, rtol=rtol, atol=0.0)


@pytest.mark.parametrize("skew", [False, True])
@pytest.mark.parametrize("m", range(1, 12))
def test_exact_cfbend_small_curvature(m, skew):
    """At small curvature, a multipole kicks as in ExactMultipole.

    The multipole of index m (quadrupole: m = 1) has a field of 1e-3 at 2 mm, in a bend
    with a radius of 1 km. Its curvilinear terms change the kick by a relative amount of
    order h x = 2e-6.
    """
    h, radius = 1.0e-3, 2.0e-3
    zero = [0.0] * (m + 1)
    k = list(zero)
    k[m] = factorial(m) * 1.0e-3 / radius**m
    kwargs = dict(ds=0.1, int_order=2, mapsteps=10)

    def change(element, k_normal_0):
        k_normal = [k_normal_0] + (zero if skew else k)[1:]
        k_skew = k if skew else zero
        with_multipole = track(
            element(k_normal=k_normal, k_skew=k_skew, **kwargs), radius
        )
        without = track(
            element(k_normal=[k_normal_0] + zero[1:], k_skew=zero, **kwargs), radius
        )
        return with_multipole - without

    # the kick: changes of px and py, relative to their largest value
    momenta = [PHASE_SPACE_COLUMNS.index(c) for c in ["momentum_x", "momentum_y"]]
    bend = change(elements.ExactCFbend, h)[:, momenta]
    straight = change(elements.ExactMultipole, 0.0)[:, momenta]
    error = np.abs(bend - straight).max(axis=0) / np.abs(straight).max(axis=0)

    rtol = 1.0e-5
    if Config.precision == "SINGLE":
        # FIXME: the bend map rounds x to multiples of eps * rc in single precision,
        # which changes the kick by up to 0.1; with #1776, the error is up to 4e-6
        # and rtol = 1e-5 holds
        # https://github.com/BLAST-ImpactX/impactx/issues/1771
        rtol = 0.3
    assert np.all(error < rtol), error
