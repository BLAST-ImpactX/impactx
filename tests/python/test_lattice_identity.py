#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""A lattice holds elements the way a Python list holds its items.

Appending stores the object given, not a copy, so the caller keeps a handle to the very
element that is tracked. Lookup and membership go by identity.
"""

import gc

import pytest

from impactx import RefPart, elements, push, reverse


def drifts(count):
    """A lattice of ``count`` named drifts, and the names in order."""
    names = [f"d{i}" for i in range(count)]
    lattice = elements.KnownElementsList(
        [elements.Drift(ds=0.1 * (i + 1), name=name) for i, name in enumerate(names)]
    )
    return lattice, names


def names_of(lattice):
    return [element.name for element in lattice]


def test_the_lattice_holds_the_object_given():
    q = elements.Quad(ds=0.3, k=2.0)
    lattice = elements.KnownElementsList()
    lattice.append(q)

    assert lattice[0] is q

    # retuning through either handle is one element
    q.k = 3.0
    assert lattice[0].k == 3.0
    lattice[0].k = 4.0
    assert q.k == 4.0


def test_the_same_element_twice_is_one_element_at_two_positions():
    q = elements.Quad(ds=0.3, k=2.0)
    lattice = elements.KnownElementsList([q, q])

    assert lattice[0] is lattice[1] is q
    q.k = 7.0
    assert [element.k for element in lattice] == [7.0, 7.0]

    # constructing (or copying) one per position gives independent elements
    independent = elements.KnownElementsList(
        [elements.Quad(ds=0.3, k=2.0), elements.Quad(ds=0.3, k=2.0)]
    )
    independent[0].k = 9.0
    assert independent[1].k == 2.0


def test_elements_and_their_python_subclass_survive_in_the_lattice():
    """The lattice keeps the element alive, including its Python subclass and attributes."""

    class Tagged(elements.Drift):
        def __init__(self, tag):
            super().__init__(ds=1.0)
            self.tag = tag

    lattice = elements.KnownElementsList()
    lattice.append(Tagged("kept"))
    lattice.append(elements.Quad(ds=0.3, k=2.0, name="only-in-lattice"))
    other = elements.KnownElementsList([lattice[1]])
    gc.collect()

    assert type(lattice[0]) is Tagged
    assert lattice[0].tag == "kept"
    assert lattice[1].name == "only-in-lattice"
    # two lattices can share one element
    assert other[0] is lattice[1]


@pytest.mark.parametrize(
    "source",
    [
        pytest.param(list, id="list"),
        pytest.param(elements.KnownElementsList, id="lattice"),
        pytest.param(lambda items: (item for item in items), id="generator"),
    ],
)
def test_constructor_takes_any_iterable(source):
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)

    lattice = elements.KnownElementsList(source([q, d]))

    assert len(lattice) == 2
    assert lattice[0] is q
    assert lattice[1] is d


def test_constructor_takes_a_single_element():
    q = elements.Quad(ds=0.3, k=2.0)

    assert elements.KnownElementsList(q)[0] is q


def test_indexing_and_iteration():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d, q])

    assert lattice[-1] is q
    assert lattice[-2] is d
    for index in (3, -4):
        with pytest.raises(IndexError):
            lattice.__getitem__(index)

    forward = list(lattice)
    assert forward[0] is q and forward[1] is d and forward[2] is q
    backward = list(reversed(lattice))
    assert backward[0] is q and backward[1] is d and backward[2] is q


def test_list_operations():
    q = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    m = elements.Marker("m")
    lattice = elements.KnownElementsList([q])

    lattice[0] = d
    assert lattice[0] is d

    lattice.insert(0, q)
    assert lattice[0] is q and lattice[1] is d
    lattice.insert(99, m)  # an out-of-range position clamps, as for a list
    assert lattice[-1] is m

    del lattice[0]
    assert lattice[0] is d
    with pytest.raises(IndexError):
        del lattice[5]

    assert lattice.pop_back() is m
    assert lattice.pop_back() is d
    assert len(lattice) == 0
    with pytest.raises(IndexError):
        lattice.pop_back()


def test_lookup_goes_by_identity():
    """Two elements with equal parameters are still two different elements."""

    q = elements.Quad(ds=0.3, k=2.0)
    twin = elements.Quad(ds=0.3, k=2.0)
    d = elements.Drift(ds=1.0)
    lattice = elements.KnownElementsList([q, d, q])

    assert q in lattice
    assert twin not in lattice
    assert "not an element" not in lattice
    assert lattice.count(q) == 2
    assert lattice.count(twin) == 0
    assert lattice.index(q) == 0
    with pytest.raises(ValueError):
        lattice.index(twin)

    # removing takes the first occurrence and leaves the element at its other position
    lattice.remove(q)
    assert list(lattice) == [d, q]
    with pytest.raises(ValueError):
        lattice.remove(twin)


def test_non_elements_are_rejected_and_change_nothing():
    q = elements.Quad(ds=0.3, k=2.0)
    lattice = elements.KnownElementsList([q])

    with pytest.raises(TypeError):
        lattice.append("not an element")
    with pytest.raises(TypeError):
        lattice.extend([elements.Drift(ds=1.0), "not an element"])
    with pytest.raises(TypeError):
        lattice[0:1] = [elements.Drift(ds=1.0), "not an element"]

    assert list(lattice) == [q]


@pytest.mark.parametrize(
    ("count", "key"),
    [
        (5, slice(1, 4)),
        (5, slice(None, None, 2)),
        (5, slice(-2, None)),
        (5, slice(None, None, -2)),
        (6, slice(4, 1, -1)),
        (3, slice(10, None)),
    ],
)
def test_slice_read_and_delete_match_a_list(count, key):
    lattice, reference = drifts(count)

    part = lattice[key]
    assert type(part) is elements.KnownElementsList
    assert names_of(part) == reference[key]
    # a slice is a new lattice over the same elements
    for element in part:
        assert element in lattice

    del lattice[key]
    del reference[key]
    assert names_of(lattice) == reference


@pytest.mark.parametrize(
    ("key", "replacements"),
    [
        (slice(0, 3), 3),
        (slice(0, 3), 5),
        (slice(1, 4), 1),
        (slice(None), 2),
        (slice(2, 2), 2),
        (slice(2, 4), 0),
        (slice(None, None, 2), 3),
    ],
)
def test_slice_assignment_matches_a_list(key, replacements):
    lattice, reference = drifts(6)
    new = [elements.Drift(ds=0.2, name=f"n{i}") for i in range(replacements)]

    lattice[key] = new
    reference[key] = [element.name for element in new]

    assert names_of(lattice) == reference
    if new:
        assert new[0] in lattice


def test_extended_slice_assignment_requires_matching_length():
    lattice, reference = drifts(4)

    with pytest.raises(ValueError):
        lattice[::2] = [elements.Drift(ds=1.0)]

    assert names_of(lattice) == reference


def test_free_functions_act_on_the_element_given():
    d = elements.Drift(ds=1.0)
    reverse(d)
    assert d.ds == -1.0

    seen = []
    hooked = elements.Programmable()
    hooked.ref_particle = lambda refpart: seen.append(refpart.s)

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(100.0)
    push(ref, hooked)

    assert seen == [0.0]
