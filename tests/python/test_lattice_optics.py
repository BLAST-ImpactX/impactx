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
    lattice.append(elements.Programmable(name="foobar"))

    # Ensure that the calculation asserts
    with pytest.raises(RuntimeError):
        lattice.transfer_map(ref)

    # Now the user explicitly assumes that undefined maps are identity maps
    R = lattice.transfer_map(ref, fallback_identity_map=True)
    assert np.allclose(R.to_numpy(), R_expected, rtol=rtol, atol=atol)


def test_programmable_linear_hook():
    """When Programmable.linear_map hook is set, it must be used by
    transfer_map. When it is unset, transfer_map must raise."""
    from impactx import Map6x6

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(100.0)

    # A user-supplied linear map: a thin quad kick with strength k_user.
    k_user = 0.7
    R_user = np.eye(6)
    R_user[1, 0] = -k_user
    R_user[3, 2] = -k_user

    def my_map(_refpart):
        out = Map6x6()
        for i in range(6):
            for j in range(6):
                out[i + 1, j + 1] = float(R_user[i, j])
        return out

    # Without hook: must raise. (extend() copies the element into the
    # lattice's variant, so any subsequent modification of the local
    # Python handle does not propagate into the lattice.)
    lattice_no_hook = elements.KnownElementsList()
    lattice_no_hook.extend(
        [elements.Drift(ds=0.1), elements.Programmable(), elements.Drift(ds=0.1)]
    )
    with pytest.raises(RuntimeError):
        lattice_no_hook.transfer_map(ref)

    # With hook set before insertion: drift(0.1) * R_user * drift(0.1).
    prog = elements.Programmable()
    prog.linear_map = my_map
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Drift(ds=0.1), prog, elements.Drift(ds=0.1)])
    R = lattice.transfer_map(ref).to_numpy()

    pt = ref.pt
    bg2 = pt * pt - 1.0
    D = np.eye(6)
    D[0, 1] = 0.1
    D[2, 3] = 0.1
    D[4, 5] = 0.1 / bg2
    R_expected = D @ R_user @ D

    tol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-10
    assert np.max(np.abs(R - R_expected)) < tol
