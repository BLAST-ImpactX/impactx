#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
"""The ``unit`` argument selects a convention and is an integer everywhere."""

import pytest

from impactx import elements

# minimal constructor arguments per element that takes an integer ``unit``
UNIT_ELEMENTS = {
    "ChrPlasmaLens": dict(ds=0.3, k=1.0),
    "ChrQuad": dict(ds=0.3, k=1.0),
    "ExactCFbend": dict(ds=0.3, k_normal=[1.0], k_skew=[0.0]),
    "ExactMultipole": dict(ds=0.3, k_normal=[1.0], k_skew=[0.0]),
    "ExactQuad": dict(ds=0.3, k=1.0),
    "QuadEdge": dict(k=1.0),
    "SoftSolenoid": dict(
        ds=1.0, bscale=1.0, cos_coefficients=[1.0], sin_coefficients=[0.0]
    ),
    "TaperedPL": dict(k=1.0, taper=0.0),
}


@pytest.mark.parametrize("name", sorted(UNIT_ELEMENTS))
@pytest.mark.parametrize("unit", [0, 1])
def test_unit_accepts_integers(name, unit):
    """Both conventions can be selected, and the value is stored as given."""

    element = getattr(elements, name)(unit=unit, **UNIT_ELEMENTS[name])

    assert element.unit == unit


@pytest.mark.parametrize("name", sorted(UNIT_ELEMENTS))
def test_unit_rejects_floats(name):
    """A float is not a convention: it is rejected rather than truncated."""

    with pytest.raises(TypeError):
        getattr(elements, name)(unit=1.7, **UNIT_ELEMENTS[name])
