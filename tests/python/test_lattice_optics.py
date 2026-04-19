#!/usr/bin/env python3
#
# Copyright 2022-2024 The ImpactX Community
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
    lattice.append(elements.TaperedPL(k=0, taper=0, unit=0))

    # Ensure that the calculation asserts
    with pytest.raises(RuntimeError):
        lattice.transfer_map(ref)

    # Now the user explicitly assumes that undefined maps are identity maps
    R = lattice.transfer_map(ref, fallback_identity_map=True)
    assert np.allclose(R.to_numpy(), R_expected, rtol=rtol, atol=atol)
