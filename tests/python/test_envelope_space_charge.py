#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import pytest

from impactx import ImpactX, distribution, elements

# beam and lattice parameters, shared by all runs below
KIN_ENERGY_MEV = 250.0
BUNCH_CHARGE_C = 1.0e-9
DS = 2.0  # drift length in m

# the beam moments that a drift under space charge changes
MOMENTS = ["sig_x", "sig_y", "sig_t", "emittance_x", "emittance_y", "emittance_t"]


def _soft_quad():
    """A soft-edge quadrupole, whose map is integrated in ``mapsteps`` steps per slice.

    Halving its slice transport re-integrates it on a finer grid, so unlike a drift it does
    not compose exactly: ``M(ds/2) M(ds/2) != M(ds)``. That makes it the element that can
    tell a split slice apart from an unsplit one.
    """
    return elements.SoftQuadrupole(
        name="sq1",
        ds=DS,
        gscale=1.0,
        cos_coefficients=[2.0],  # a flat on-axis profile of 1
        sin_coefficients=[0.0],
        mapsteps=8,
        nslice=4,
    )


def _envelope_sim(space_charge, intensity=BUNCH_CHARGE_C, lattice=None):
    """A lattice under the given space charge model, set up for envelope tracking."""
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = space_charge
    sim.slice_step_diagnostics = False
    sim.diagnostics = False

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)

    distr = distribution.Kurth6D(
        lambdaX=4.472135955e-4,
        lambdaY=4.472135955e-4,
        lambdaT=9.12241869e-7,
        lambdaPx=1.0e-4,
        lambdaPy=1.0e-4,
        lambdaPt=0.0,
    )
    sim.init_envelope(ref, distr, intensity)

    if lattice is None:
        lattice = [elements.Drift(name="d1", ds=DS, nslice=10)]
    sim.lattice.extend(lattice)

    return sim


@pytest.mark.parametrize("space_charge", ["Gauss3D", "Gauss2p5D", "2p5D"])
def test_particle_only_space_charge_models_are_rejected(space_charge):
    """Envelope tracking implements the 2D and 3D models only.

    The models it does not implement have to be rejected rather than tracked without any
    space charge at all, which would report a space-charge-free answer for a run the user
    asked to include it.
    """
    sim = _envelope_sim(space_charge)

    try:
        with pytest.raises(RuntimeError, match="only supported with particle tracking"):
            sim.track_envelope()
    finally:
        sim.finalize()


def _run(space_charge, intensity):
    """Track the soft quadrupole and return its beam moments."""
    sim = _envelope_sim(space_charge, intensity=intensity, lattice=[_soft_quad()])

    try:
        sim.track_envelope()
        return sim.envelope.beam_moments(sim.beam.ref)
    finally:
        sim.finalize()


def test_zero_intensity_is_the_same_as_no_space_charge():
    """A zero-charge beam takes no collective kick at all.

    The kick is linear in the intensity, so at zero it leaves the covariance matrix
    untouched and the run has to reproduce a space-charge-free one exactly, rather than
    merely to within the error of a split element transport.
    """
    zero_intensity = _run("3D", intensity=0.0)
    no_space_charge = _run("false", intensity=0.0)

    for key in MOMENTS:
        assert zero_intensity[key] == no_space_charge[key], (
            f"{key}: {zero_intensity[key]} != {no_space_charge[key]} at zero intensity"
        )
