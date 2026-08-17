#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import impactx
from impactx import ImpactX, distribution


def waterbag(lambda_pt=2.0e-3):
    """A 2 GeV electron beam, optionally without any energy spread."""
    return distribution.Waterbag(
        lambdaX=3.9984884770e-5,
        lambdaY=3.9984884770e-5,
        lambdaT=1.0e-3,
        lambdaPx=2.6623538760e-5,
        lambdaPy=2.6623538760e-5,
        lambdaPt=lambda_pt,
        muxpx=-0.846574929020762,
        muypy=0.846574929020762,
        mutpt=0.0,
    )


def test_beam_moments_single_particle():
    """
    A single particle has vanishing rms emittances, for which the Courant-Snyder
    (Twiss) functions are undefined.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(1.0e-9, waterbag(), 1)

    moments = sim.beam.beam_moments()

    # the beam is a point in phase space: no extent, no emittance
    for key in ("sigma_x", "sigma_y", "sigma_t", "emittance_x", "emittance_y"):
        assert moments[key] == 0.0

    # ... so its Twiss functions do not exist, in any plane
    for key in ("alpha_x", "alpha_y", "alpha_t", "beta_x", "beta_y", "beta_t"):
        assert math.isnan(moments[key])

    # its position is still perfectly well defined
    assert math.isfinite(moments["mean_x"])
    assert moments["mean_x"] == moments["min_x"] == moments["max_x"]

    sim.finalize()


def test_beam_moments_without_energy_spread():
    """
    A beam without energy spread has no longitudinal emittance, and hence no
    longitudinal Twiss functions. The transverse planes stay unaffected.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(1.0e-9, waterbag(lambda_pt=0.0), 1000)

    moments = sim.beam.beam_moments()

    assert moments["sigma_pt"] == 0.0
    assert moments["emittance_t"] == 0.0
    assert math.isnan(moments["beta_t"])
    assert math.isnan(moments["alpha_t"])

    assert moments["sigma_t"] > 0.0
    assert moments["emittance_x"] > 0.0
    assert math.isfinite(moments["beta_x"])
    assert math.isfinite(moments["alpha_x"])

    sim.finalize()


def test_envelope_beam_moments_without_energy_spread():
    """
    The covariance-matrix (envelope) moments treat a vanishing emittance the same
    way as the particle moments do.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    envelope = impactx.create_envelope(waterbag(lambda_pt=0.0), 1.0e-9)
    moments = envelope.beam_moments(ref)

    assert moments["emittance_t"] == 0.0
    assert math.isnan(moments["beta_t"])
    assert math.isnan(moments["alpha_t"])

    assert moments["emittance_x"] > 0.0
    assert math.isfinite(moments["beta_x"])
    assert math.isfinite(moments["alpha_x"])

    sim.finalize()
