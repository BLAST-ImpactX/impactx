#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, RefPart, elements


def test_lattice_linear_map():
    """Calculate the linear transfer map of a lattice."""

    # Create reference particle
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    # Create a valid test lattice (all elements define a linear transfer map)
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Sbend(name="bend1", ds=1.0, rc=10.0, nslice=2),
            elements.Drift(name="drift1", ds=2.0, nslice=2),
            elements.Quad(name="quad1", ds=0.5, k=1.0, nslice=2),
            elements.Drift(name="drift2", ds=1.0, nslice=2),
        ]
    )

    # Expected result (matrix multiplication)
    R_expected = np.array(
        [
            [3.74670546e-01, 2.54005827e00, 0, 0, 0, -2.34864805e-01],
            [-4.76219076e-01, -5.59489407e-01, 0, 0, 0, 3.20646253e-02],
            [0, 0, 1.64872127e00, 6.59488508e00, 0, 0],
            [0, 0, 5.21095305e-01, 2.69091188e00, 0, 0],
            [9.98334297e-02, 4.99583537e-02, 0, 0, 1.00000000e00, -1.66466013e-03],
            [0, 0, 0, 0, 0, 1.00000000e00],
        ]
    )

    if Config.precision == "SINGLE":
        atol = 1.0e-7
        rtol = 5.0e-5
    else:
        atol = 0.0
        rtol = 1.0e-8

    # Calculate Linear Transfer Map
    R = lattice.transfer_map(ref)
    assert np.allclose(R.to_numpy(), R_expected, rtol=rtol, atol=atol)

    # Check unexpected/unsupported options
    with pytest.raises(RuntimeError):
        lattice.transfer_map(ref, order="invalid")

    # Create a lattice with an element that does not define a linear transfer map
    lattice.append(elements.Programmable(name="foobar"))

    # Ensure that the calculation asserts
    with pytest.raises(RuntimeError):
        lattice.transfer_map(ref)

    # Now the user explicitly assumes that undefined maps are identity maps
    R = lattice.transfer_map(ref, fallback_identity_map=True)
    assert np.allclose(R.to_numpy(), R_expected, rtol=rtol, atol=atol)


def _single_element_lattices():
    """One-slice elements whose linear map depends on the reference push."""
    return [
        # plain control: map is independent of the reference push
        elements.Drift(ds=0.3),
        # map is evaluated at the post-advance reference momentum
        elements.ShortRF(V=5.0e6, freq=1.3e9, phase=-90.0),
        elements.ChrAcc(ds=0.3, ez=1.0e6, bz=1.0),
        # map is integrated during the reference push
        elements.RFCavity(
            ds=0.3,
            escale=1.0,
            freq=1.3e9,
            phase=-90.0,
            cos_coefficients=[2.0],
            sin_coefficients=[0.0],
        ),
        elements.SoftQuadrupole(
            ds=0.3, gscale=1.0, cos_coefficients=[2.0], sin_coefficients=[0.0]
        ),
        elements.SoftSolenoid(
            ds=0.3, bscale=1.0, cos_coefficients=[2.0], sin_coefficients=[0.0]
        ),
    ]


@pytest.mark.parametrize("element", _single_element_lattices())
def test_element_transfer_map_matches_lattice(element):
    """An element's own transfer map must equal the map of a lattice that
    contains only that element (both are the single-slice map)."""

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    if Config.precision == "SINGLE":
        rtol, atol = 5.0e-5, 1.0e-7
    else:
        rtol, atol = 1.0e-8, 0.0

    lattice = elements.KnownElementsList()
    lattice.append(element)

    R_element = element.transfer_map(ref).to_numpy()
    R_lattice = lattice.transfer_map(ref).to_numpy()

    assert np.allclose(R_element, R_lattice, rtol=rtol, atol=atol)
    # a linear map is never singular; the map-integrating elements used to
    # return an all-zero matrix before their reference push had ever run
    assert not np.allclose(R_element, 0.0)


def test_element_transfer_map_keeps_reference_particle():
    """Evaluating an element's transfer map must not advance the caller's
    reference particle."""

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    s_before, pt_before = ref.s, ref.pt

    elements.ShortRF(V=5.0e6, freq=1.3e9, phase=-90.0).transfer_map(ref)
    elements.Drift(ds=0.3).transfer_map(ref)

    assert ref.s == s_before
    assert ref.pt == pt_before
