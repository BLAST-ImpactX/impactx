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
NPART = 1000
DS = 2.0  # drift length in m
NSLICE = 3
PERIODS = 4

# the tracking state every model has to agree on: the last period that was tracked, and
# one step per slice of every element in every period
LAST_PERIOD = PERIODS - 1
TOTAL_STEPS = PERIODS * NSLICE


def _distribution():
    return distribution.Kurth6D(
        lambdaX=4.472135955e-4,
        lambdaY=4.472135955e-4,
        lambdaT=9.12241869e-7,
        lambdaPx=1.0e-4,
        lambdaPy=1.0e-4,
        lambdaPt=0.0,
    )


def _sim():
    """A single drift, tracked for several periods, without any collective effects."""
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = "false"
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.periods = PERIODS

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)

    sim.lattice.extend([elements.Drift(name="d1", ds=DS, nslice=NSLICE)])

    return sim


def _track(model):
    """Run one tracking model and return its tracking state afterwards."""
    sim = _sim()

    try:
        if model == "particles":
            sim.add_particles(BUNCH_CHARGE_C, _distribution(), NPART)
            sim.track_particles()
        elif model == "envelope":
            sim.init_envelope(sim.beam.ref, _distribution(), BUNCH_CHARGE_C)
            sim.track_envelope()
        else:
            sim.track_reference(sim.beam.ref)

        return sim.tracking_period, sim.tracking_step
    finally:
        sim.finalize()


@pytest.mark.parametrize("model", ["particles", "envelope", "reference"])
def test_tracking_state_after_a_multi_period_run(model):
    """Every tracking model traverses the lattice the same way.

    The period is the last one that was tracked, not one past it, and the step counter
    has advanced once per slice of every element in every period. These are the values
    a hook and ``reduced_beam_characteristics_final`` report, so the three models have to
    agree on them.
    """
    period, step = _track(model)

    assert period == LAST_PERIOD, (
        f"{model}: tracked period {period}, expected {LAST_PERIOD}"
    )
    assert step == TOTAL_STEPS, f"{model}: {step} steps, expected {TOTAL_STEPS}"
