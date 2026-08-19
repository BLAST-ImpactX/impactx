#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import numpy as np
import pytest

import amrex.space3d as amr
from impactx import ImpactX, Map6x6, elements

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


def _new_sim(strang_split):
    """A simulation with space charge enabled and the deterministic beam loaded."""
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

    return sim


def _drift_map(sim, ds):
    """The 6x6 transport map of a drift of length ``ds``, as a LinearMap would apply it."""
    betgam2 = sim.beam.ref.pt**2 - 1.0

    R = Map6x6.identity()
    R[1, 2] = ds
    R[3, 4] = ds
    R[5, 6] = ds / betgam2

    return R


def _run(strang_split, nslice):
    """Track the beam through a single drift under space charge.

    :param strang_split: second-order Strang split (True) or first-order composition
    :param nslice: number of slices through the drift
    :return: the beam characteristics after the drift
    """
    sim = _new_sim(strang_split)

    sim.lattice.extend([elements.Drift(name="d1", ds=DS, nslice=nslice)])

    sim.track_particles()

    moments = sim.beam.beam_moments()
    sim.finalize()

    return moments


def _chain(sim, sliceable, n_elements):
    """A chain of ``n_elements`` drift-equivalent elements adding up to a length of ``DS``.

    With ``sliceable``, plain drifts are used, which are subdivided into half transports.
    Otherwise the same drift is expressed as a ``LinearMap``, which has a finite ``ds`` but
    cannot be sliced, so the kick is halved instead. Refining the integration then means
    shorter elements rather than more slices.
    """
    ds = DS / n_elements

    if sliceable:
        return [elements.Drift(ds=ds) for _ in range(n_elements)]

    R = _drift_map(sim, ds)
    return [elements.LinearMap(R=R, ds=ds) for _ in range(n_elements)]


def _run_chain(strang_split, sliceable, n_elements):
    """Track the beam through a chain of elements under space charge."""
    sim = _new_sim(strang_split)

    sim.lattice.extend(_chain(sim, sliceable, n_elements))

    sim.track_particles()

    moments = sim.beam.beam_moments()
    sim.finalize()

    return moments


def _roundtrip(strang_split, sliceable, n_elements=4):
    """Track the chain forward, then back through the same elements reversed.

    :return: the largest relative deviation of the returned beam from the initial one
    """
    sim = _new_sim(strang_split)

    forward = _chain(sim, sliceable, n_elements)
    backward = _chain(sim, sliceable, n_elements)
    for element in reversed(backward):
        element.reverse()

    initial = sim.beam.beam_moments()

    sim.lattice.extend(forward + backward)
    sim.track_particles()

    final = sim.beam.beam_moments()
    sim.finalize()

    return max(
        abs(final[key] / initial[key] - 1.0)
        for key in ["sig_x", "sig_y", "sig_px", "sig_py", "emittance_x", "emittance_y"]
    )


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


def _observed_order_chain(strang_split, sliceable, n_elements=4):
    """Estimate the order of convergence in the element length, @see _observed_order."""
    coarse = _run_chain(strang_split, sliceable, n_elements)["sig_x"]
    medium = _run_chain(strang_split, sliceable, 2 * n_elements)["sig_x"]
    fine = _run_chain(strang_split, sliceable, 4 * n_elements)["sig_x"]

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


def test_element_that_cannot_be_sliced_is_second_order():
    """An element that cannot be sliced halves the kick, which is second order as well."""
    order = _observed_order_chain(strang_split=True, sliceable=False)

    assert order > 1.7, f"halved kick converges at order {order:.2f}, expected ~2"


def test_element_that_cannot_be_sliced_is_first_order_when_disabled():
    """Disabling the split falls back to first order for those elements, too."""
    order = _observed_order_chain(strang_split=False, sliceable=False)

    assert order < 1.3, (
        f"first-order composition converges at order {order:.2f}, expected ~1"
    )


@pytest.mark.parametrize(
    "sliceable", [True, False], ids=["halved_transport", "halved_kick"]
)
def test_split_is_time_symmetric(sliceable):
    """Tracking on through the reversed elements restores the initial beam.

    Both splits are their own time-reverse, so the collective kicks of the return pass cancel
    those of the forward pass exactly, down to round-off.
    """
    residual = _roundtrip(strang_split=True, sliceable=sliceable)

    assert residual < 1.0e-10, f"beam did not return: {residual:.3e}"


@pytest.mark.parametrize(
    "sliceable", [True, False], ids=["transport_first", "kick_first"]
)
def test_first_order_composition_is_not_time_symmetric(sliceable):
    """The first-order composition does not return the beam this way.

    It inverts exactly when the lattice is traversed backwards, since the tracking loop then
    mirrors the order of kick and transport. Reversing the elements and tracking on forward
    keeps that order, which leaves the transport between the two kicks.
    """
    residual = _roundtrip(strang_split=False, sliceable=sliceable)

    assert residual > 1.0e-4, (
        f"expected a residual from the first-order composition: {residual:.3e}"
    )
