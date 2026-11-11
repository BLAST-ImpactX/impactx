#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""The linear map of an RFCavity is continuous in its length and does not depend on
the position of the cavity in the beamline.

The integration evaluates the on-axis field profile, which lies on 0 <= z <= ds, at
positions accumulated step by step from s - sedge. Round-off can put an evaluation at
an end of the profile a few ulp outside of it.
"""

import numpy as np
import pytest

from impactx import Config, RefPart, elements

DS = 1.31879807

# an on-axis field profile with a finite slope at both ends of the cavity
PROFILE = dict(
    escale=62.0,
    freq=1.3e9,
    phase=85.5,
    cos_coefficients=[2.0, 0.3, -0.1],
    sin_coefficients=[0.0, 0.1, -0.03],
)

if Config.precision == "SINGLE":
    dtype = np.float32
    rtol = 1.0e-4
else:
    dtype = np.float64
    rtol = 1.0e-11


def electron_1GeV():
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    return ref


def assert_maps_equal(maps):
    """All maps agree with the first, entry by entry to rtol, and entries close to
    zero to a tenth of rtol of the largest entry"""
    maps = np.array(maps, dtype=np.float64)
    np.testing.assert_allclose(
        maps,
        np.broadcast_to(maps[0], maps.shape),
        rtol=rtol,
        atol=0.1 * rtol * np.abs(maps[0]).max(),
    )


@pytest.mark.parametrize("mapsteps", [3, 10, 64])
def test_map_is_continuous_in_ds(mapsteps):
    """Changing ds by a few ulp does not change the map"""
    ds = dtype(DS)
    for _ in range(8):
        ds = np.nextafter(ds, dtype(0.0))
    maps = []
    for _ in range(17):
        el = elements.RFCavity(ds=float(ds), mapsteps=mapsteps, **PROFILE)
        maps.append(el.transfer_map(electron_1GeV()).to_numpy().copy())
        ds = np.nextafter(ds, dtype(10.0))
    assert_maps_equal(maps)


@pytest.mark.parametrize("nslice", [2, 3, 4, 5, 8])
def test_map_is_independent_of_s(nslice):
    """The map of a cavity is the same at any position s in the beamline

    The first slice starts at s - sedge = 0 exactly, so at least two slices.
    """
    lattice = elements.KnownElementsList()
    lattice.append(elements.RFCavity(ds=DS, nslice=nslice, **PROFILE))
    maps = []
    for s0 in [0.0, 0.1, 3.7, 12.345, 123.456]:
        ref = electron_1GeV()
        ref.s = s0
        maps.append(lattice.transfer_map(ref).to_numpy().copy())
    assert_maps_equal(maps)
