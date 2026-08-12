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
NSLICE = 4
QM_EEV = -1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV


def _deterministic_beam():
    """A fixed, reproducible beam in s-coordinates relative to the reference particle.

    Sampled with a fixed seed instead of ``ImpactX.add_particles``, because the AMReX RNG
    stream advances between runs in the same process. Every run in this test must start
    from bit-identical particles, otherwise the sampling masks the tested effect.
    """
    rng = np.random.default_rng(seed=42)

    lambda_x = 4.472135955e-4
    lambda_t = 9.12241869e-7

    x = rng.normal(0.0, lambda_x, NPART)
    y = rng.normal(0.0, lambda_x, NPART)
    t = rng.normal(0.0, lambda_t, NPART)
    # a cold beam: the expansion below is driven by space charge alone
    px = np.zeros(NPART)
    py = np.zeros(NPART)
    pt = np.zeros(NPART)

    return x, y, t, px, py, pt


def _drift_beam_particles(pge, pti, refpart):
    """Push beam particles as a drift, equivalent to elements.Drift.

    Mirrors examples/fodo_programmable/run_fodo_programmable.py. ``nslice`` is read live
    on every call, which is the documented Programmable contract: "for an element with
    nslice > 1, the pushes and maps refer to a single ds/nslice slice".
    """
    soa = pti.soa().to_xp()

    x = soa.real["position_x"]
    y = soa.real["position_y"]
    t = soa.real["position_t"]
    px = soa.real["momentum_x"]
    py = soa.real["momentum_y"]
    pt = soa.real["momentum_t"]

    slice_ds = pge.ds / pge.nslice

    betgam2 = refpart.pt**2 - 1.0

    x[:] += slice_ds * px[:]
    y[:] += slice_ds * py[:]
    t[:] += (slice_ds / betgam2) * pt[:]


def _drift_ref_particle(pge, refpart):
    """Push the reference particle as a drift."""
    slice_ds = pge.ds / pge.nslice
    step = slice_ds / (refpart.pt**2 - 1.0) ** 0.5

    refpart.x = refpart.x + step * refpart.px
    refpart.y = refpart.y + step * refpart.py
    refpart.z = refpart.z + step * refpart.pz
    refpart.t = refpart.t - step * refpart.pt
    refpart.s = refpart.s + slice_ds


def _run(programmable, space_charge, nslice=NSLICE):
    """Track a cold beam through a single drift and return its beam characteristics.

    :param programmable: use a Programmable element instead of elements.Drift
    :param space_charge: enable the 3D space-charge solver
    :param nslice: number of slices through the element
    """
    sim = ImpactX()

    sim.n_cell = [16, 16, 16]
    sim.particle_shape = 2
    sim.space_charge = "3D" if space_charge else False
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

    if programmable:
        element = elements.Programmable(name="d1")
        element.ds = DS
        element.nslice = nslice
        element.beam_particles = lambda pti, refpart: _drift_beam_particles(
            element, pti, refpart
        )
        element.ref_particle = lambda refpart: _drift_ref_particle(element, refpart)
    else:
        element = elements.Drift(name="d1", ds=DS, nslice=nslice)

    sim.lattice.extend([element])

    sim.track_particles()

    rbc = sim.beam.reduced_beam_characteristics()
    sim.finalize()

    return rbc


def test_programmable_receives_collective_kick():
    """A Programmable element with a finite ``ds`` must receive the collective kick.

    A Programmable element carries its own ``ds`` and ``nslice`` instead of deriving from
    the ``Thick`` mixin. It is therefore not eligible for the Strang split, but it is
    eligible for the kick itself: it must not be tracked as pure optics.
    """
    prog_sc = _run(programmable=True, space_charge=True)
    prog_no_sc = _run(programmable=True, space_charge=False)

    # this beam starts cold, so any transverse momentum can only come from the kick
    assert prog_no_sc["sig_px"] == 0.0
    assert prog_sc["sig_px"] > 0.0, (
        "the Programmable element received no collective kick"
    )

    # and the kick has to visibly expand the beam
    sc_effect = prog_sc["sig_x"] / prog_no_sc["sig_x"] - 1.0
    assert sc_effect > 0.05, f"space charge expanded sig_x by only {sc_effect:.2%}"


def test_programmable_transports_like_an_equivalent_drift():
    """A Programmable drift must transport the same length as elements.Drift, kick or not.

    A Programmable element cannot be subdivided by overriding its slice count. Its push is a
    Python callback that reads ``ds`` and ``nslice`` from the user's own object, while the
    lattice holds a copy, so the override is invisible to it and each half-map would transport
    a full slice. It therefore uses the first-order composition, which agrees with the
    Strang-split Drift only in the converged limit.
    """
    # without collective effects both take the identical, unsplit path
    drift = _run(programmable=False, space_charge=False)
    prog = _run(programmable=True, space_charge=False)
    for key in ["sig_x", "sig_y", "sig_t", "emittance_x", "emittance_y"]:
        assert math.isclose(drift[key], prog[key], rel_tol=1e-12), (
            f"{key}: Drift={drift[key]:.8e} vs Programmable={prog[key]:.8e}"
        )

    # with space charge the two integrators differ at O(1/nslice), but must converge to the
    # same physics. A double-length transport would not shrink with nslice.
    coarse = _relative_difference(NSLICE)
    fine = _relative_difference(4 * NSLICE)
    assert fine < coarse, (
        f"not converging: {coarse:.4%} at nslice={NSLICE} -> {fine:.4%}"
    )
    assert fine < 0.02, (
        f"Drift and Programmable differ by {fine:.4%} at nslice={4 * NSLICE}"
    )


def _relative_difference(nslice):
    """Relative sig_x difference between Drift and Programmable under space charge."""
    drift = _run(programmable=False, space_charge=True, nslice=nslice)
    prog = _run(programmable=True, space_charge=True, nslice=nslice)
    return abs(prog["sig_x"] / drift["sig_x"] - 1.0)
