#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import numpy as np

import amrex.space3d as amr
from impactx import ImpactX, elements

# beam and lattice parameters, shared by all runs below
KIN_ENERGY_MEV = 250.0
BUNCH_CHARGE_C = 1.0e-9
NPART = 10000
DS = 2.0  # drift length in m
QM_EEV = -1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV


def _deterministic_beam():
    """A fixed, reproducible beam in s-coordinates relative to the reference particle.

    Sampled with a fixed seed instead of ``ImpactX.add_particles``, because the AMReX RNG
    stream advances between runs in the same process. Every run in this test must start
    from bit-identical particles, otherwise the sampling masks the tested effect.

    The beam is diverging: with a finite ``px``, the leading half transport of the Strang
    split moves the particles before the first collective kick is applied.
    """
    rng = np.random.default_rng(seed=42)

    lambda_x = 4.472135955e-4
    lambda_px = 1.0e-4
    lambda_t = 9.12241869e-7

    x = rng.normal(0.0, lambda_x, NPART)
    y = rng.normal(0.0, lambda_x, NPART)
    t = rng.normal(0.0, lambda_t, NPART)
    px = rng.normal(0.0, lambda_px, NPART)
    py = rng.normal(0.0, lambda_px, NPART)
    pt = np.zeros(NPART)

    return x, y, t, px, py, pt


def _run(strang_split, nslice):
    """Track the beam through a single drift under space charge.

    :param strang_split: second-order Strang split (True) or first-order composition
    :param nslice: number of slices through the drift
    :return: the beam characteristics after the drift
    """
    sim = ImpactX()

    sim.n_cell = [16, 16, 16]
    sim.particle_shape = 2
    sim.space_charge = "3D"
    sim.strang_split = strang_split
    sim.dynamic_size = True
    sim.prob_relative = [3.0]
    sim.slice_step_diagnostics = False
    sim.diagnostics = False

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)
    ref.z = 0.0

    if amr.ParallelDescriptor.IOProcessor():
        x, y, t, px, py, pt = _deterministic_beam()
        sim.beam.add_n_particles(
            x, y, t, px, py, pt, QM_EEV, bunch_charge=BUNCH_CHARGE_C
        )

    sim.lattice.extend([elements.Drift(name="d1", ds=DS, nslice=nslice)])

    sim.track_particles()

    moments = sim.beam.beam_moments()
    sim.finalize()

    return moments


def _observed_order(strang_split, nslice=4):
    """Estimate the order of convergence in the slice length.

    Runs the same beam at ``nslice``, ``2 * nslice`` and ``4 * nslice`` and compares the two
    successive differences. Halving the slice length shrinks the error by ``2**order``, so
    the ratio of the differences gives the order without needing an exact solution.
    """
    coarse = _run(strang_split, nslice)["sig_x"]
    medium = _run(strang_split, 2 * nslice)["sig_x"]
    fine = _run(strang_split, 4 * nslice)["sig_x"]

    return math.log2(abs(coarse - medium) / abs(medium - fine))


def _relative_difference(nslice):
    """Relative sig_x difference between the second-order and first-order composition."""
    split = _run(strang_split=True, nslice=nslice)["sig_x"]
    first = _run(strang_split=False, nslice=nslice)["sig_x"]

    return abs(split / first - 1.0)


def test_strang_split_is_second_order():
    """The default composition converges with the square of the slice length."""
    order = _observed_order(strang_split=True)

    assert order > 1.7, f"Strang split converges at order {order:.2f}, expected ~2"


def test_first_order_composition():
    """Disabling the split falls back to first-order convergence."""
    order = _observed_order(strang_split=False)

    assert order < 1.3, (
        f"first-order composition converges at order {order:.2f}, expected ~1"
    )


def test_both_compositions_converge_to_the_same_result():
    """The two compositions differ at a given slicing, but not in the converged limit."""
    coarse = _relative_difference(nslice=8)
    fine = _relative_difference(nslice=32)

    assert coarse > 1.0e-3, "algo.strang_split had no effect on the tracked beam"

    # the gap is dominated by the first-order error, so it shrinks with the slice length
    assert fine < coarse / 2.5, f"not converging: {coarse:.4%} -> {fine:.4%}"
    assert fine < 3.0e-3, f"the two compositions differ by {fine:.4%} at nslice=32"
