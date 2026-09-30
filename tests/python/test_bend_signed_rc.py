#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""``signed_rc(ref)``: the signed radius of curvature of the reference orbit in a bend."""

import math

import pytest

from impactx import ImpactX, RefPart, elements


def electron(kin_energy_MeV=2.0e3):
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)
    return ref


def bends(rc):
    """One bend of each kind with a reference orbit of radius rc, bending the same way"""
    ds = 0.5
    phi_deg = math.degrees(ds / rc)
    return [
        elements.Sbend(ds=ds, rc=rc),
        elements.ExactSbend(ds=ds, phi=phi_deg),
        elements.CFbend(ds=ds, rc=rc, k=0.1),
        elements.ExactCFbend(ds=ds, k_normal=[1.0 / rc, 0.1], k_skew=[0.0, 0.0]),
    ]


@pytest.mark.parametrize("rc", [10.0, -10.0])
def test_all_bends_agree(rc):
    for bend in bends(rc):
        assert bend.signed_rc(electron()) == pytest.approx(rc, rel=1e-12), bend


def test_radius_from_a_field_scales_with_the_rigidity():
    B = 0.5  # T
    for kin_energy_MeV in [1.0e3, 2.0e3]:
        ref = electron(kin_energy_MeV)
        rc = ref.rigidity_Tm / B

        # the bend angle only sets the magnitude of the angle, the field the radius
        exact_sbend = elements.ExactSbend(ds=0.5, phi=10.0, B=B)
        assert exact_sbend.signed_rc(ref) == pytest.approx(rc, rel=1e-12)

        # unit=1: the dipole coefficient is a field in T
        exact_cfbend = elements.ExactCFbend(ds=0.5, k_normal=[B], k_skew=[0.0], unit=1)
        assert exact_cfbend.signed_rc(ref) == pytest.approx(rc, rel=1e-12)


@pytest.mark.parametrize(
    "element",
    [
        elements.Drift(ds=0.5),
        elements.Quad(ds=0.5, k=1.0),
        elements.ExactMultipole(ds=0.5, k_normal=[0.1], k_skew=[0.0]),
        elements.ThinDipole(theta=1.0, rc=10.0),
        elements.DipEdge(psi=0.1, rc=10.0, g=0.05),
    ],
)
def test_only_bends_with_synchrotron_radiation_have_it(element):
    assert not hasattr(element, "signed_rc")


@pytest.mark.parametrize("rc", [10.0, -10.0])
def test_the_center_of_curvature_is_at_negative_signed_rc(rc):
    """A bend with a positive signed_rc turns the reference orbit towards negative x"""
    for bend in bends(rc):
        sim = ImpactX()
        sim.particle_shape = 2
        sim.slice_step_diagnostics = False
        sim.diagnostics = False
        sim.init_grids()
        ref = sim.beam.ref
        ref.set_species("electron").set_kin_energy_MeV(2.0e3)

        signed_rc = bend.signed_rc(ref)
        sim.lattice.append(bend)
        sim.track_reference(ref)

        # the orbit stays on the circle around (x, z) = (-signed_rc, 0)
        assert math.copysign(1.0, ref.x) == -math.copysign(1.0, signed_rc), bend
        assert math.hypot(ref.x + signed_rc, ref.z) == pytest.approx(
            abs(signed_rc), rel=1e-12
        ), bend

        sim.finalize()
