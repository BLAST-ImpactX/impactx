#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Test that zero-length (thin) lattice elements apply a transverse beam pipe
aperture, just like the finite-length (thick) elements do.
"""

import numpy as np
import pytest

from impactx import ImpactX, distribution, elements

# beam sizes, chosen so that a 1 sigma aperture cuts a sizable fraction
LAMBDA_X = 1.0e-3
LAMBDA_Y = 2.0e-3

# half-apertures used for the "particles are lost" checks
APERTURE_X = LAMBDA_X
APERTURE_Y = LAMBDA_Y

NPART = 10000


def thin_elements(aperture_x, aperture_y):
    """One instance of each thin element, all with zero strength.

    Zero strength keeps the transverse coordinates untouched by the push, so
    the surviving particles must lie inside the aperture ellipse exactly.
    """
    return {
        "Buncher": elements.Buncher(
            V=0.0, k=0.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "DipEdge": elements.DipEdge(
            psi=0.0, rc=1.0, g=0.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "Kicker": elements.Kicker(
            xkick=0.0, ykick=0.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "Multipole": elements.Multipole(
            multipole=1,
            K_normal=0.0,
            K_skew=0.0,
            aperture_x=aperture_x,
            aperture_y=aperture_y,
        ),
        "NonlinearLens": elements.NonlinearLens(
            knll=0.0, cnll=1.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "PRot": elements.PRot(
            phi_in=0.0, phi_out=0.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "PlaneXYRot": elements.PlaneXYRot(
            angle=0.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "QuadEdge": elements.QuadEdge(
            k=0.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "ShortRF": elements.ShortRF(
            V=0.0, freq=1.0e6, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "TaperedPL": elements.TaperedPL(
            k=0.0, taper=0.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
        "ThinDipole": elements.ThinDipole(
            theta=0.0, rc=1.0, aperture_x=aperture_x, aperture_y=aperture_y
        ),
    }


THIN_ELEMENT_NAMES = sorted(thin_elements(0.0, 0.0).keys())


def _track(element):
    """Track a Gaussian beam through a single element and return its container."""
    sim = ImpactX()

    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    distr = distribution.Gaussian(
        lambdaX=LAMBDA_X,
        lambdaY=LAMBDA_Y,
        lambdaT=1.0e-3,
        lambdaPx=0.0,
        lambdaPy=0.0,
        lambdaPt=0.0,
    )
    sim.add_particles(1.0e-9, distr, NPART)

    pc = sim.beam
    assert pc.total_number_of_particles() == NPART

    sim.lattice.append(element)
    sim.track_particles()

    return sim, pc


@pytest.mark.parametrize("element_name", THIN_ELEMENT_NAMES)
def test_thin_element_aperture_scrapes(element_name):
    """A tight aperture on a thin element removes the particles outside it."""
    element = thin_elements(APERTURE_X, APERTURE_Y)[element_name]
    sim, pc = _track(element)

    n_final = pc.total_number_of_particles()
    assert n_final < NPART, f"{element_name}: no particle was scraped"

    # the surviving particles lie inside the aperture ellipse
    df = pc.to_df()
    r2 = (df["position_x"] / APERTURE_X) ** 2 + (df["position_y"] / APERTURE_Y) ** 2
    assert np.all(r2 <= 1.0 + 1.0e-5), (
        f"{element_name}: a particle outside the aperture survived"
    )

    sim.finalize()


@pytest.mark.parametrize("element_name", THIN_ELEMENT_NAMES)
def test_thin_element_aperture_off_by_default(element_name):
    """An aperture of zero (the default) does not scrape any particle."""
    element = thin_elements(0.0, 0.0)[element_name]
    sim, pc = _track(element)

    assert pc.total_number_of_particles() == NPART, (
        f"{element_name}: particles were lost without an aperture"
    )

    sim.finalize()
