#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

NPART = 10000
SIGMA = 1.0e-3  # transverse rms beam size in m
HALF_APERTURE = 1.0e-3  # in m


def _track(aperture):
    """
    Track a round Gaussian beam through a single aperture element.

    :param aperture: the thin collimator to track through
    :return: number of surviving particles, their max |x| and max |y| in m
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("proton").set_kin_energy_MeV(250.0)

    distr = distribution.Gaussian(
        lambdaX=SIGMA,
        lambdaY=SIGMA,
        lambdaT=1.0e-3,
        lambdaPx=1.0e-4,
        lambdaPy=1.0e-4,
        lambdaPt=1.0e-3,
    )
    sim.add_particles(1.0e-9, distr, NPART)

    sim.lattice.append(aperture)
    sim.track_particles()

    pc = sim.particle_container()
    n_surviving = pc.total_number_of_particles()

    max_x, max_y = 0.0, 0.0
    if n_surviving > 0:
        xmin, ymin, _, xmax, ymax, _ = pc.min_and_max_positions()
        max_x = max(abs(xmin), abs(xmax))
        max_y = max(abs(ymin), abs(ymax))

    sim.finalize()

    return n_surviving, max_x, max_y


def test_aperture_disabled_transmits_everything():
    """
    Test that an Aperture with both half-apertures disabled is a no-op for the
    default "transmit" action, for both the constructor and the setters.
    """
    n_ctor, _, _ = _track(elements.Aperture(aperture_x=0.0, aperture_y=0.0))
    assert n_ctor == NPART

    collimator = elements.Aperture(aperture_x=HALF_APERTURE, aperture_y=HALF_APERTURE)
    collimator.aperture_x = 0.0
    collimator.aperture_y = 0.0
    n_setter, _, _ = _track(collimator)
    assert n_setter == NPART


def test_aperture_disabled_plane_is_a_slit():
    """
    Test that disabling a single plane keeps the boundary of the other plane:
    a rectangular aperture with only aperture_y set is a horizontal slit.
    """
    n_slit, max_x, max_y = _track(
        elements.Aperture(aperture_x=0.0, aperture_y=HALF_APERTURE)
    )
    n_both, max_x_both, _ = _track(
        elements.Aperture(aperture_x=HALF_APERTURE, aperture_y=HALF_APERTURE)
    )

    # the vertical boundary is enforced in both cases; the boundary itself is
    # transmitting, so allow for the roundoff of the inverted-aperture scaling
    on_edge = HALF_APERTURE * (1.0 + 1.0e-9)
    assert 0 < n_slit < NPART
    assert max_y <= on_edge

    # the disabled horizontal plane does not restrict x, so the slit keeps
    # everything the fully bounded aperture keeps, and more
    assert n_slit > n_both
    assert max_x > HALF_APERTURE
    assert max_x_both <= on_edge


def test_aperture_disabled_plane_absorbs_unbounded():
    """
    Test the documented corner of the "zero disables" convention: for the
    "absorb" action a disabled plane makes the absorbing domain unbounded, so
    disabling both planes absorbs the whole beam.
    """
    n_surviving, _, _ = _track(
        elements.Aperture(aperture_x=0.0, aperture_y=0.0, action="absorb")
    )
    assert n_surviving == 0
