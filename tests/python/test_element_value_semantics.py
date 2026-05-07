#!/usr/bin/env python3
#
# Copyright 2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Tests for value-based __eq__, __hash__, and isclose() on lattice elements.
"""

import inspect

import pytest

import amrex.space3d as amr
from impactx import elements

# ----------------------------------------------------------------------
# Basic equality and hash
# ----------------------------------------------------------------------


def test_eq_same_params():
    a = elements.Drift(ds=1.0, name="d1")
    b = elements.Drift(ds=1.0, name="d1")
    assert a == b
    assert hash(a) == hash(b)


def test_eq_different_params():
    a = elements.Drift(ds=1.0, name="d1")
    b = elements.Drift(ds=2.0, name="d1")
    assert a != b


def test_eq_different_types():
    # Drift and ExactDrift have the same constructor kwargs but are
    # distinct element types; they must not compare equal.
    a = elements.Drift(ds=1.0)
    b = elements.ExactDrift(ds=1.0)
    assert a != b


def test_eq_foreign_type_returns_notimplemented():
    # Direct API check: __eq__ returns NotImplemented for foreign operands.
    a = elements.Drift(ds=1.0)
    assert a.__eq__("foo") is NotImplemented
    # Python's reflected fallback should make `==` evaluate to False.
    assert (a == "foo") is False
    assert a != "foo"


def test_set_membership():
    s = {elements.Drift(ds=1.0, name="d1")}
    assert elements.Drift(ds=1.0, name="d1") in s
    assert elements.Drift(ds=2.0, name="d1") not in s


def test_dict_key_usage():
    a = elements.Drift(ds=1.0, name="d1")
    b = elements.Drift(ds=1.0, name="d1")
    d = {a: "first"}
    assert d[b] == "first"


# ----------------------------------------------------------------------
# isclose: scalars
# ----------------------------------------------------------------------


def test_isclose_accepts_small_perturbation():
    a = elements.Drift(ds=1.0, name="d1")
    b = elements.Drift(ds=1.0 + 1e-15, name="d1")
    assert a.isclose(b)
    assert not (a == b)


def test_isclose_rejects_large_perturbation():
    a = elements.Drift(ds=1.0)
    b = elements.Drift(ds=1.001)
    assert not a.isclose(b)


def test_isclose_custom_rtol():
    a = elements.Drift(ds=1.0)
    b = elements.Drift(ds=1.0 + 1e-9)
    # default rtol=1e-12 rejects this
    assert not a.isclose(b)
    # loosened rtol accepts it
    assert a.isclose(b, rtol=1e-6)


def test_isclose_with_atol_only():
    a = elements.Drift(ds=0.0)
    b = elements.Drift(ds=1e-10)
    # rtol=0, atol=1e-9 should accept
    assert a.isclose(b, rtol=0.0, atol=1e-9)
    # rtol=0, atol=0 should reject (math.isclose default behavior)
    assert not a.isclose(b, rtol=0.0, atol=0.0)


def test_isclose_different_types_returns_false():
    a = elements.Drift(ds=1.0)
    b = elements.ExactDrift(ds=1.0)
    assert not a.isclose(b)


def test_isclose_foreign_operand():
    a = elements.Drift(ds=1.0)
    assert not a.isclose("foo")


# ----------------------------------------------------------------------
# isclose: list-of-floats (Fourier coefficients, polygon vertices)
# ----------------------------------------------------------------------


def _make_rfcavity(cos_coefs, sin_coefs):
    return elements.RFCavity(
        ds=0.5,
        escale=1e6,
        freq=1.3e9,
        phase=-90.0,
        cos_coefficients=cos_coefs,
        sin_coefficients=sin_coefs,
        nslice=2,
        mapsteps=10,
        name="rf",
    )


def test_eq_list_of_floats():
    a = _make_rfcavity([1.0, 2.0, 3.0], [0.0, 0.1, 0.2])
    b = _make_rfcavity([1.0, 2.0, 3.0], [0.0, 0.1, 0.2])
    assert a == b
    assert hash(a) == hash(b)


def test_isclose_list_of_floats_perturbed():
    a = _make_rfcavity([1.0, 2.0, 3.0], [0.0, 0.1, 0.2])
    b = _make_rfcavity([1.0, 2.0, 3.0 + 1e-15], [0.0, 0.1, 0.2])
    assert a.isclose(b)
    assert not (a == b)


def test_isclose_list_of_floats_rejected():
    a = _make_rfcavity([1.0, 2.0, 3.0], [0.0, 0.1, 0.2])
    b = _make_rfcavity([1.0, 2.0, 3.001], [0.0, 0.1, 0.2])
    assert not a.isclose(b)


def test_isclose_list_length_mismatch():
    a = _make_rfcavity([1.0, 2.0], [0.0, 0.1])
    b = _make_rfcavity([1.0, 2.0, 3.0], [0.0, 0.1, 0.2])
    assert not a.isclose(b)
    assert a != b


# ----------------------------------------------------------------------
# SmallMatrix-valued entries (LinearMap)
# ----------------------------------------------------------------------


def test_eq_linearmap_same_matrix():
    R = amr.SmallMatrix_6x6_F_SI1_double.identity()
    a = elements.LinearMap(R=R, ds=0.5, name="lm")
    b = elements.LinearMap(R=R, ds=0.5, name="lm")
    assert a == b
    assert hash(a) == hash(b)


def test_isclose_linearmap_perturbed_matrix():
    R1 = amr.SmallMatrix_6x6_F_SI1_double.identity()
    R2 = amr.SmallMatrix_6x6_F_SI1_double.identity()
    # SmallMatrix_..._SI1_... uses 1-based indexing.
    R2[1, 1] = 1.0 + 1e-15
    a = elements.LinearMap(R=R1, ds=0.5, name="lm")
    b = elements.LinearMap(R=R2, ds=0.5, name="lm")
    assert a.isclose(b)
    assert not (a == b)


# ----------------------------------------------------------------------
# NaN handling
# ----------------------------------------------------------------------


def test_nan_not_equal():
    a = elements.Drift(ds=float("nan"))
    b = elements.Drift(ds=float("nan"))
    assert a != b
    assert not a.isclose(b)


def test_nan_hashes_consistently():
    # Hashing must not raise on NaN values.
    a = elements.Drift(ds=float("nan"))
    hash(a)


# ----------------------------------------------------------------------
# Sweep over every patched element class
# ----------------------------------------------------------------------


def _patched_element_classes():
    for name in dir(elements):
        if name.startswith("_"):
            continue
        cls = getattr(elements, name)
        if not inspect.isclass(cls):
            continue
        if any(
            cls is element_list_type
            for element_list_type in elements._ELEMENT_LIST_TYPES
        ):
            continue
        if not hasattr(cls, "to_dict"):
            continue
        yield name, cls


def test_every_element_class_is_patched():
    found_any = False
    for name, cls in _patched_element_classes():
        found_any = True
        # __eq__ replaced
        assert cls.__eq__ is not object.__eq__, f"{name}: __eq__ not patched"
        # __hash__ kept (not None) so instances stay hashable
        assert cls.__hash__ is not None, f"{name}: __hash__ is None"
        # isclose attached
        assert callable(getattr(cls, "isclose", None)), f"{name}: isclose missing"
    assert found_any, "no element classes were discovered"


def test_marker_self_equality():
    # Marker is a thin element with minimal state — sanity check that the
    # patched methods work for a stripped-down class.
    a = elements.Marker(name="m")
    b = elements.Marker(name="m")
    assert a == b
    assert hash(a) == hash(b)
    assert a.isclose(b)


@pytest.mark.parametrize(
    "factory",
    [
        lambda: elements.Drift(ds=1.0, name="d"),
        lambda: elements.Marker(name="m"),
        lambda: elements.Quad(ds=0.3, k=2.5, name="q"),
        lambda: elements.Sbend(ds=0.5, rc=5.0, name="s"),
    ],
    ids=["Drift", "Marker", "Quad", "Sbend"],
)
def test_self_equality_and_hash(factory):
    a = factory()
    b = factory()
    assert a == b
    assert hash(a) == hash(b)
    assert a.isclose(b)


# ----------------------------------------------------------------------
# Container-level value semantics: KnownElementsList & FilteredElementsList
# ----------------------------------------------------------------------


def _build_lattice():
    lat = elements.KnownElementsList()
    lat.append(elements.Drift(ds=1.0, name="d1"))
    lat.append(elements.Quad(ds=0.3, k=2.5, name="q1"))
    lat.append(elements.Drift(ds=2.0, name="d2"))
    return lat


def test_lattice_eq_same_contents():
    a = _build_lattice()
    b = _build_lattice()
    assert a == b


def test_lattice_eq_different_length():
    a = _build_lattice()
    b = _build_lattice()
    b.append(elements.Marker(name="m"))
    assert a != b


def test_lattice_eq_different_element():
    a = _build_lattice()
    b = elements.KnownElementsList()
    b.append(elements.Drift(ds=1.0, name="d1"))
    b.append(elements.Quad(ds=2.5, k=2.5, name="q1"))  # ds differs
    b.append(elements.Drift(ds=2.0, name="d2"))
    assert a != b


def test_lattice_eq_with_python_list():
    # Duck-typed cross-type equality with a plain Python list of elements.
    a = _build_lattice()
    py_list = [
        elements.Drift(ds=1.0, name="d1"),
        elements.Quad(ds=0.3, k=2.5, name="q1"),
        elements.Drift(ds=2.0, name="d2"),
    ]
    assert a == py_list
    assert py_list == a  # symmetric


def test_lattice_eq_with_filtered_view():
    # Cross-type duck-typed equality between a KnownElementsList and a
    # FilteredElementsList covering the same indices.
    a = _build_lattice()
    drifts = a.select(kind="Drift")
    expected = elements.KnownElementsList()
    expected.append(elements.Drift(ds=1.0, name="d1"))
    expected.append(elements.Drift(ds=2.0, name="d2"))
    assert drifts == expected
    assert expected == drifts  # symmetric across container types


def test_lattice_eq_filtered_vs_filtered():
    a = _build_lattice()
    b = _build_lattice()
    assert a.select(kind="Drift") == b.select(kind="Drift")
    assert a.select(kind="Drift") != b.select(kind="Quad")


def test_lattice_eq_foreign_type():
    a = _build_lattice()
    assert a != "lattice"
    with pytest.raises(ValueError):
        assert a != 42
    with pytest.raises(ValueError):
        assert a != None  # noqa: E711


def test_lattice_isclose_with_perturbation():
    a = _build_lattice()
    b = elements.KnownElementsList()
    b.append(elements.Drift(ds=1.0 + 1e-15, name="d1"))
    b.append(elements.Quad(ds=0.3 + 1e-15, k=2.5, name="q1"))
    b.append(elements.Drift(ds=2.0, name="d2"))
    assert a.isclose(b)
    assert not (a == b)


def test_lattice_isclose_rejects_large_perturbation():
    a = _build_lattice()
    b = elements.KnownElementsList()
    b.append(elements.Drift(ds=1.001, name="d1"))
    b.append(elements.Quad(ds=0.3, k=2.5, name="q1"))
    b.append(elements.Drift(ds=2.0, name="d2"))
    assert not a.isclose(b)
    assert a.isclose(b, rtol=1e-2)


def test_lattice_isclose_length_mismatch():
    a = _build_lattice()
    b = _build_lattice()
    b.append(elements.Marker(name="m"))
    assert not a.isclose(b)


def test_lattice_isclose_foreign_type():
    a = _build_lattice()
    assert not a.isclose("lattice")
    with pytest.raises(ValueError):
        assert not a.isclose(None)


def test_lattice_isclose_filtered_view():
    a = _build_lattice()
    drifts = a.select(kind="Drift")
    expected = [
        elements.Drift(ds=1.0 + 1e-15, name="d1"),
        elements.Drift(ds=2.0, name="d2"),
    ]
    assert drifts.isclose(expected)


# ----------------------------------------------------------------------
# isclose ignore_attributes
# ----------------------------------------------------------------------


def test_isclose_ignore_name():
    a = elements.Drift(ds=1.0, name="alpha")
    b = elements.Drift(ds=1.0, name="beta")
    assert not a.isclose(b)
    assert a.isclose(b, ignore_attributes=["name"])


def test_isclose_ignore_single_string():
    # Convenience: a single string is treated as a one-element set.
    a = elements.Drift(ds=1.0, name="alpha")
    b = elements.Drift(ds=1.0, name="beta")
    assert a.isclose(b, ignore_attributes="name")


def test_isclose_ignore_type_across_variants():
    # Drift and ExactDrift share the same parameter set; ignoring "type"
    # lets them compare on common fields.
    a = elements.Drift(
        ds=1.0,
        dx=0.001,
        dy=0.002,
        rotation=0.05,
        aperture_x=0.01,
        aperture_y=0.02,
        name="d",
    )
    b = elements.ExactDrift(
        ds=1.0,
        dx=0.001,
        dy=0.002,
        rotation=0.05,
        aperture_x=0.01,
        aperture_y=0.02,
        name="d",
    )
    assert not a.isclose(b)
    assert a.isclose(b, ignore_attributes=["type"])


def test_isclose_ignore_type_keys_still_must_match():
    # Drift and Quad have different keys (Quad has "k") even after
    # dropping "type"; isclose must still return False.
    a = elements.Drift(ds=1.0, name="x")
    b = elements.Quad(ds=1.0, k=2.5, name="x")
    assert not a.isclose(b, ignore_attributes=["type"])


def test_isclose_ignore_multiple():
    a = elements.Drift(ds=1.0, name="alpha", dx=0.001)
    b = elements.Drift(ds=1.0, name="beta", dx=0.999)
    assert a.isclose(b, ignore_attributes={"name", "dx"})


def test_isclose_ignore_attributes_empty():
    # Empty iterable behaves like None (no filtering).
    a = elements.Drift(ds=1.0, name="alpha")
    b = elements.Drift(ds=1.0, name="beta")
    assert not a.isclose(b, ignore_attributes=[])


def test_lattice_isclose_ignore_attributes_forwarded():
    # ignore_attributes should be forwarded to each element's isclose.
    a = elements.KnownElementsList()
    a.append(elements.Drift(ds=1.0, name="d_a"))
    a.append(elements.Quad(ds=0.3, k=2.5, name="q_a"))

    b = elements.KnownElementsList()
    b.append(elements.Drift(ds=1.0, name="d_b"))
    b.append(elements.Quad(ds=0.3, k=2.5, name="q_b"))

    assert not a.isclose(b)
    assert a.isclose(b, ignore_attributes=["name"])


def test_lattice_isclose_ignore_type_across_variants():
    # Mixed-variant lattices: ignoring "type" lets a linear lattice
    # compare equal to its exact-variant counterpart on common fields.
    a = elements.KnownElementsList()
    a.append(elements.Drift(ds=1.0, name="d1"))
    a.append(elements.Drift(ds=2.0, name="d2"))

    b = elements.KnownElementsList()
    b.append(elements.ExactDrift(ds=1.0, name="d1"))
    b.append(elements.ExactDrift(ds=2.0, name="d2"))

    assert not a.isclose(b)
    assert a.isclose(b, ignore_attributes=["type"])


def test_lattice_uses_identity_hash():
    # Containers do not get a value-based __hash__: the inherited identity
    # hash is preserved (so internal weak-reference tracking of filtered
    # selections keeps working). Two value-equal containers therefore may
    # hash differently — they're not meant as dict/set keys for value-based
    # deduplication.
    a = _build_lattice()
    b = _build_lattice()
    assert a == b
    assert hash(a) != hash(b)  # identity-based, distinct instances
    # Filtered views are also hashable (identity-based).
    drifts_a = a.select(kind="Drift")
    drifts_b = a.select(kind="Drift")
    assert drifts_a == drifts_b
    assert hash(drifts_a) != hash(drifts_b)
