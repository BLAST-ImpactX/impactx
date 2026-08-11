#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
Test enhanced KnownElementsList functionality with a select method.

This module tests the select method that supports:
- Filtering by element type (kind)
- Filtering by element name (including regex patterns)
- OR-based filtering with lists/tuples
- AND-based chaining of filters
- Reference preservation for modification
"""

import pytest

from impactx import elements
from impactx.extensions.KnownElementsList import (
    FilteredElementsList as _FilteredElementsList,
)


def test_filtered_elements_list_exported_on_elements_module():
    assert elements.FilteredElementsList is _FilteredElementsList
    assert elements.FilteredElementsList.__module__ == elements.__name__


def test_basic_indexing():
    """Test basic integer indexing functionality."""

    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
        ]
    )

    # Test integer indexing
    assert type(lattice[0]) is elements.Drift
    assert lattice[0].name == "drift1"
    assert type(lattice[1]) is elements.Quad
    assert lattice[1].name == "quad1"
    assert type(lattice[3]) is elements.Sbend
    assert lattice[3].name == "bend1"

    # Test modification through integer indexing
    lattice[0].ds = 5.0
    assert lattice[0].ds == 5.0  # Original element should be modified

    # Test modification through integer indexing (use ds which is modifiable)
    lattice[3].ds = 15.0
    assert lattice[3].ds == 15.0  # Original element should be modified


def test_single_value_filtering():
    """Test filtering by single kind and name values."""

    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
            elements.ChrQuad(name="quad2", ds=0.3, k=-1.0),
        ]
    )

    # Test filtering by kind (string-based)
    drift_elements = lattice.select(kind="Drift")
    assert len(drift_elements) == 2
    assert all(type(el) is elements.Drift for el in drift_elements)
    assert drift_elements[0].name == "drift1"
    assert drift_elements[1].name == "drift2"

    quad_elements = lattice.select(kind=r".*Quad")
    assert len(quad_elements) == 2
    assert all(type(el) in [elements.Quad, elements.ChrQuad] for el in quad_elements)

    bend_elements = lattice.select(kind="Sbend")
    assert len(bend_elements) == 1
    assert bend_elements[0].name == "bend1"

    # Test filtering by kind (type-based)
    drift_elements_type = lattice.select(kind=elements.Drift)
    assert len(drift_elements_type) == 2
    assert all(type(el) is elements.Drift for el in drift_elements_type)
    assert drift_elements_type[0].name == "drift1"
    assert drift_elements_type[1].name == "drift2"

    quad_elements_type = lattice.select(kind=[elements.Quad, elements.ChrQuad])
    assert len(quad_elements_type) == 2
    assert all(
        type(el) in [elements.Quad, elements.ChrQuad] for el in quad_elements_type
    )

    bend_elements_type = lattice.select(kind=elements.Sbend)
    assert len(bend_elements_type) == 1
    assert bend_elements_type[0].name == "bend1"

    # Verify string and type filtering give same results
    assert len(drift_elements) == len(drift_elements_type)
    assert len(quad_elements) == len(quad_elements_type)
    assert len(bend_elements) == len(bend_elements_type)

    # Test filtering by name
    drift1_elements = lattice.select(name="drift1")
    assert len(drift1_elements) == 1
    assert type(drift1_elements[0]) is elements.Drift
    assert drift1_elements[0].name == "drift1"

    quad1_elements = lattice.select(name="quad1")
    assert len(quad1_elements) == 1
    assert type(quad1_elements[0]) is elements.Quad

    # Test regex filtering by name
    drift_regex_elements = lattice.select(name=r"drift\d+")
    assert len(drift_regex_elements) == 2
    assert all(type(el) is elements.Drift for el in drift_regex_elements)
    assert drift_regex_elements[0].name == "drift1"
    assert drift_regex_elements[1].name == "drift2"

    # Test manipulation of regex-selected drift elements
    drift_regex_elements[0].ds = 10.0
    drift_regex_elements[1].ds = 20.0
    assert lattice[0].ds == 10.0  # Original element should be modified
    assert lattice[2].ds == 20.0  # Original element should be modified

    quad_regex_elements = lattice.select(name=r"quad\d+")
    assert len(quad_regex_elements) == 2
    assert all(
        type(el) in [elements.Quad, elements.ChrQuad] for el in quad_regex_elements
    )
    assert quad_regex_elements[0].name == "quad1"
    assert quad_regex_elements[1].name == "quad2"

    # Test manipulation of regex-selected quad elements
    quad_regex_elements[0].k = 5.0
    quad_regex_elements[1].k = -5.0
    assert lattice[1].k == 5.0  # Original element should be modified
    assert lattice[4].k == -5.0  # Original element should be modified

    # Test regex that matches only one element
    bend_regex_elements = lattice.select(name=r"bend\d+")
    assert len(bend_regex_elements) == 1
    assert type(bend_regex_elements[0]) is elements.Sbend
    assert bend_regex_elements[0].name == "bend1"

    # Test manipulation of regex-selected bend element
    bend_regex_elements[0].ds = 5.0
    assert lattice[3].ds == 5.0  # Original element should be modified

    # Test regex that matches no elements
    empty_regex_elements = lattice.select(name=r"nonexistent\d+")
    assert len(empty_regex_elements) == 0

    # Test regex filtering by kind with r".*Quad" pattern
    quad_kind_regex_elements = lattice.select(kind=r".*Quad")
    assert len(quad_kind_regex_elements) == 2
    assert all(
        type(el) in [elements.Quad, elements.ChrQuad] for el in quad_kind_regex_elements
    )
    assert quad_kind_regex_elements[0].name == "quad1"
    assert quad_kind_regex_elements[1].name == "quad2"

    # Test FilteredElementsList helper methods
    assert quad_kind_regex_elements.get_kinds() == [elements.ChrQuad, elements.Quad]
    assert quad_kind_regex_elements.count_by_kind("Quad") == 1
    assert quad_kind_regex_elements.count_by_kind("ChrQuad") == 1
    assert quad_kind_regex_elements.count_by_kind(elements.Drift) == 0
    assert quad_kind_regex_elements.has_kind("Quad")
    assert quad_kind_regex_elements.has_kind(elements.ChrQuad)
    assert not quad_kind_regex_elements.has_kind("Drift")

    # Test manipulation of regex-selected quad elements by kind
    quad_kind_regex_elements[0].k = 7.0
    quad_kind_regex_elements[1].k = -7.0
    assert lattice[1].k == 7.0  # Original element should be modified
    assert lattice[4].k == -7.0  # Original element should be modified

    # Test combined filtering (OR logic between different filters)
    quads_combined = lattice.select(kind=r".*Quad", name="quad1")
    assert (
        len(quads_combined) == 2
    )  # All Quad elements (2) OR element named "quad1" (1) = 2 total
    assert all(type(el) in [elements.Quad, elements.ChrQuad] for el in quads_combined)


def test_or_based_filtering():
    """Test OR-based filtering with lists and tuples."""

    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
            elements.Quad(name="quad2", ds=0.3, k=-1.0),
            elements.Marker(name="marker1"),
        ]
    )

    # Test OR filtering by multiple kinds (string list)
    drift_quad_elements = lattice.select(kind=["Drift", "Quad"])
    assert len(drift_quad_elements) == 4  # 2 drifts + 2 quads
    assert all(
        type(el) in [elements.Drift, elements.Quad] for el in drift_quad_elements
    )

    # Test OR filtering by multiple kinds (type list)
    drift_quad_type_elements = lattice.select(kind=[elements.Drift, elements.Quad])
    assert len(drift_quad_type_elements) == 4  # 2 drifts + 2 quads
    assert all(
        type(el) in [elements.Drift, elements.Quad] for el in drift_quad_type_elements
    )

    # Test OR filtering by multiple kinds (mixed string and type)
    mixed_kind_elements = lattice.select(kind=["Drift", elements.Quad])
    assert len(mixed_kind_elements) == 4  # 2 drifts + 2 quads
    assert all(
        type(el) in [elements.Drift, elements.Quad] for el in mixed_kind_elements
    )

    # Test FilteredElementsList helper methods on mixed filtering
    mixed_kinds_list = mixed_kind_elements.get_kinds()
    assert elements.Drift in mixed_kinds_list
    assert elements.Quad in mixed_kinds_list
    assert len(mixed_kinds_list) == 2
    assert mixed_kind_elements.count_by_kind("Drift") == 2
    assert mixed_kind_elements.count_by_kind("Quad") == 2
    assert mixed_kind_elements.has_kind("Drift")
    assert mixed_kind_elements.has_kind("Quad")
    assert not mixed_kind_elements.has_kind("Sbend")

    # Test OR filtering by multiple kinds (tuple)
    drift_quad_tuple = lattice.select(kind=("Drift", "Quad"))
    assert len(drift_quad_tuple) == 4
    assert all(type(el) in [elements.Drift, elements.Quad] for el in drift_quad_tuple)

    # Verify all approaches give same results
    assert len(drift_quad_elements) == len(drift_quad_type_elements)
    assert len(drift_quad_elements) == len(mixed_kind_elements)
    assert len(drift_quad_elements) == len(drift_quad_tuple)

    # Test OR filtering by multiple names
    named_elements = lattice.select(name=["drift1", "quad1", "bend1"])
    assert len(named_elements) == 3
    assert all(el.name in ["drift1", "quad1", "bend1"] for el in named_elements)

    # Test combined OR filtering (OR logic between different filters)
    combined_elements = lattice.select(kind=["Drift", "Quad"], name=["drift1", "quad1"])
    assert (
        len(combined_elements) == 4
    )  # All Drift (2) OR all Quad (2) OR drift1 (1) OR quad1 (1) = 4 total
    assert all(
        el.name in ["drift1", "drift2", "quad1", "quad2"] for el in combined_elements
    )
    assert all(type(el).__name__ in ["Drift", "Quad"] for el in combined_elements)

    # Test empty list
    empty_elements = lattice.select(kind=[])
    assert len(empty_elements) == 0

    # Test single item list (should work same as string)
    single_list_elements = lattice.select(kind=["Drift"])
    single_string_elements = lattice.select(kind="Drift")
    assert len(single_list_elements) == len(single_string_elements)
    assert all(type(el) is elements.Drift for el in single_list_elements)


def test_and_based_filtering():
    """Test AND chaining of filters with regex and modifications."""

    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
            elements.Quad(name="quad2", ds=0.3, k=-1.0),
            elements.Marker(name="marker1"),
            elements.Drift(name="drift3", ds=3.0),
            elements.Quad(name="quad3", ds=0.4, k=2.0),
        ]
    )

    # Test chaining kind filter then name filter
    drift_then_name = lattice.select(kind="Drift").select(name="drift1")
    assert len(drift_then_name) == 1
    assert type(drift_then_name[0]) is elements.Drift
    assert drift_then_name[0].name == "drift1"

    # Test chaining name filter then kind filter
    name_then_kind = lattice.select(name="quad1").select(kind="Quad")
    assert len(name_then_kind) == 1
    assert type(name_then_kind[0]) is elements.Quad
    assert name_then_kind[0].name == "quad1"

    # Test chaining with OR filters
    multi_kind_then_name = lattice.select(kind=["Drift", "Quad"]).select(name="drift1")
    assert len(multi_kind_then_name) == 1
    assert type(multi_kind_then_name[0]) is elements.Drift
    assert multi_kind_then_name[0].name == "drift1"

    # Test regex-based chaining: kind then regex name
    drift_regex = lattice.select(kind=["Drift", "Quad"]).select(name=r"drift\d+")
    assert len(drift_regex) == 3  # drift1, drift2, drift3
    assert all(type(el) is elements.Drift for el in drift_regex)
    assert all(el.name in ["drift1", "drift2", "drift3"] for el in drift_regex)

    # Test regex-based chaining: regex name then kind
    regex_drift = lattice.select(name=r"drift\d+").select(kind="Drift")
    assert len(regex_drift) == 3  # drift1, drift2, drift3
    assert all(type(el) is elements.Drift for el in regex_drift)

    # Test regex-based chaining: regex name then specific name
    specific_drift = lattice.select(name=r"drift\d+").select(name="drift1")
    assert len(specific_drift) == 1
    assert type(specific_drift[0]) is elements.Drift
    assert specific_drift[0].name == "drift1"

    # Test regex-based chaining: quad regex
    quad_regex = lattice.select(kind="Quad").select(name=r"quad\d+")
    assert len(quad_regex) == 3  # quad1, quad2, quad3
    assert all(type(el) is elements.Quad for el in quad_regex)

    # Test FilteredElementsList helper methods on chained filtering
    assert quad_regex.get_kinds() == [elements.Quad]
    assert quad_regex.count_by_kind("Quad") == 3
    assert quad_regex.count_by_kind("Drift") == 0
    assert quad_regex.has_kind("Quad")
    assert not quad_regex.has_kind("Drift")

    # Test manipulation of chained regex results
    drift_regex[0].ds = 10.0  # Modify drift1
    drift_regex[1].ds = 20.0  # Modify drift2
    drift_regex[2].ds = 30.0  # Modify drift3
    assert lattice[0].ds == 10.0  # drift1 should be modified
    assert lattice[2].ds == 20.0  # drift2 should be modified
    assert lattice[6].ds == 30.0  # drift3 should be modified

    # Test manipulation of quad regex results
    quad_regex[0].k = 5.0  # Modify quad1
    quad_regex[1].k = -5.0  # Modify quad2
    quad_regex[2].k = 10.0  # Modify quad3
    assert lattice[1].k == 5.0  # quad1 should be modified
    assert lattice[4].k == -5.0  # quad2 should be modified
    assert lattice[7].k == 10.0  # quad3 should be modified


def test_error_handling():
    """Test error handling for invalid operations."""

    lattice = elements.KnownElementsList()
    lattice.append(elements.Drift(name="drift1", ds=1.0))

    # Test invalid index
    with pytest.raises(IndexError):
        _ = lattice[10]

    # Test invalid key type
    with pytest.raises(TypeError, match="incompatible function arguments"):
        _ = lattice[1.5]

    # Test invalid kind parameter type
    with pytest.raises(TypeError, match="'kind' parameter must be a string or type"):
        _ = lattice.select(kind=123)

    # Test invalid name parameter type
    with pytest.raises(
        TypeError,
        match="'name' parameter must be a string or list/tuple of strings",
    ):
        _ = lattice.select(name=123)

    # Test invalid kind parameter with mixed types in list
    with pytest.raises(
        TypeError,
        match="'kind' parameter must be a string/type or list of strings/types",
    ):
        _ = lattice.select(kind=["Drift", 123])

    # Test invalid name parameter with mixed types in list
    with pytest.raises(TypeError, match="'name' parameter must contain strings"):
        _ = lattice.select(name=["drift1", 123])


def test_empty_operations():
    """Test operations on empty lists."""

    empty_lattice = elements.KnownElementsList()

    # Test filtering on empty list
    empty_drifts = empty_lattice.select(kind="Drift")
    assert len(empty_drifts) == 0

    empty_named = empty_lattice.select(name="test")
    assert len(empty_named) == 0

    # Test chaining on empty list
    empty_chain = empty_lattice.select(kind="Drift").select(name="test")
    assert len(empty_chain) == 0


def test_helper_methods():
    """Test the helper methods for element selection and analysis."""

    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
            elements.Quad(name="quad2", ds=0.3, k=-1.0),
        ]
    )

    # Test select method for kind filtering (replaces get_by_kind)
    drift_elements = lattice.select(kind="Drift")
    assert len(drift_elements) == 2
    assert all(type(el) is elements.Drift for el in drift_elements)

    # Test select method for name filtering (replaces get_by_name)
    drift1_elements = lattice.select(name="drift1")
    assert len(drift1_elements) == 1
    assert drift1_elements[0].name == "drift1"

    # Test get_kinds method
    kinds = lattice.get_kinds()
    assert elements.Drift in kinds
    assert elements.Quad in kinds
    assert elements.Sbend in kinds
    assert len(kinds) == 3

    # Test count_by_kind method
    assert lattice.count_by_kind("Drift") == 2
    assert lattice.count_by_kind(elements.Quad) == 2
    assert lattice.count_by_kind("Sbend") == 1
    assert lattice.count_by_kind("Marker") == 0

    # Test has_kind method
    assert lattice.has_kind("Drift")
    assert lattice.has_kind("Quad")
    assert not lattice.has_kind(elements.Marker)


def test_complex_scenarios():
    """Test complex real-world scenarios."""

    # Create a realistic accelerator lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="qf1", ds=0.3, k=1.0),
            elements.Drift(name="drift2", ds=0.5),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
            elements.Drift(name="drift3", ds=0.5),
            elements.Quad(name="qd1", ds=0.3, k=-1.0),
            elements.Drift(name="drift4", ds=1.0),
            elements.Quad(name="qf2", ds=0.3, k=1.0),
            elements.Drift(name="drift5", ds=0.5),
            elements.Sbend(name="bend2", ds=1.0, rc=10.0),
        ]
    )

    # Scenario 1: Find all focusing quadrupoles (string-based OR logic between filters)
    focusing_quads = lattice.select(kind="Quad", name=["qf1", "qf2"])
    assert (
        len(focusing_quads) == 3
    )  # All Quad elements (qf1, qd1, qf2) OR elements named qf1/qf2
    assert all(el.name in ["qf1", "qd1", "qf2"] for el in focusing_quads)

    # Scenario 1b: Find all focusing quadrupoles (type-based OR logic between filters)
    focusing_quads_type = lattice.select(kind=elements.Quad, name=["qf1", "qf2"])
    assert (
        len(focusing_quads_type) == 3
    )  # All Quad elements (qf1, qd1, qf2) OR elements named qf1/qf2
    assert all(el.name in ["qf1", "qd1", "qf2"] for el in focusing_quads_type)

    # Scenario 1c: Find all focusing quadrupoles (mixed string/type OR logic)
    focusing_quads_mixed = lattice.select(
        kind=["Quad", elements.Sbend], name=["qf1", "qf2"]
    )
    assert (
        len(focusing_quads_mixed) == 5
    )  # All Quad/Sbend elements OR elements named qf1/qf2
    assert all(
        el.name in ["qf1", "qd1", "qf2", "bend1", "bend2"]
        for el in focusing_quads_mixed
    )

    # Scenario 2: Find all defocusing quadrupoles (string-based OR logic between filters)
    defocusing_quads = lattice.select(kind="Quad", name="qd1")
    assert len(defocusing_quads) == 3  # All Quad elements OR elements named qd1
    assert all(el.name in ["qf1", "qd1", "qf2"] for el in defocusing_quads)

    # Scenario 2b: Find all defocusing quadrupoles (type-based OR logic between filters)
    defocusing_quads_type = lattice.select(kind=elements.Quad, name="qd1")
    assert len(defocusing_quads_type) == 3  # All Quad elements OR elements named qd1
    assert all(el.name in ["qf1", "qd1", "qf2"] for el in defocusing_quads_type)

    # Verify string and type approaches give same results
    assert len(focusing_quads) == len(focusing_quads_type)
    assert len(defocusing_quads) == len(defocusing_quads_type)

    # Scenario 3: Find all quadrupoles and modify their strength (string-based)
    all_quads = lattice.select(kind="Quad")
    assert len(all_quads) == 3

    # Modify all quadrupoles
    for quad in all_quads:
        quad.k *= 1.1  # Increase strength by 10%

    # Verify modifications affected original elements
    assert lattice[1].k == pytest.approx(1.1)  # qf1
    assert lattice[5].k == pytest.approx(-1.1)  # qd1
    assert lattice[7].k == pytest.approx(1.1)  # qf2

    # Scenario 3b: Find all quadrupoles and modify their strength (type-based)
    all_quads_type = lattice.select(kind=elements.Quad)
    assert len(all_quads_type) == 3

    # Modify all quadrupoles using type-based filtering
    for quad in all_quads_type:
        quad.k *= 1.2  # Increase strength by 20%

    # Verify modifications affected original elements
    assert lattice[1].k == pytest.approx(1.32)  # qf1 (1.1 * 1.2)
    assert lattice[5].k == pytest.approx(-1.32)  # qd1 (-1.1 * 1.2)
    assert lattice[7].k == pytest.approx(1.32)  # qf2 (1.1 * 1.2)

    # Scenario 4: Find all bends and modify their length (string-based)
    all_bends = lattice.select(kind="Sbend")
    assert len(all_bends) == 2

    for bend in all_bends:
        bend.ds *= 0.9  # Decrease length by 10%

    # Verify modifications
    assert lattice[3].ds == pytest.approx(0.9)  # bend1
    assert lattice[9].ds == pytest.approx(0.9)  # bend2

    # Scenario 4b: Find all bends and modify their length (type-based)
    all_bends_type = lattice.select(kind=elements.Sbend)
    assert len(all_bends_type) == 2

    for bend in all_bends_type:
        bend.ds *= 0.8  # Decrease length by 20%

    # Verify modifications
    assert lattice[3].ds == pytest.approx(0.72)  # bend1 (0.9 * 0.8)
    assert lattice[9].ds == pytest.approx(0.72)  # bend2 (0.9 * 0.8)

    # Scenario 5: Chain filters to find specific elements (string-based)
    first_focusing_quad = lattice.select(kind="Quad").select(name="qf1")[0]
    assert first_focusing_quad.name == "qf1"

    # Scenario 5b: Chain filters to find specific elements (type-based)
    first_focusing_quad_type = lattice.select(kind=elements.Quad).select(name="qf1")[0]
    assert first_focusing_quad_type.name == "qf1"

    # Scenario 5c: Chain filters with mixed string/type
    first_focusing_quad_mixed = lattice.select(kind=elements.Quad).select(name="qf1")[0]
    assert first_focusing_quad_mixed.name == "qf1"
    assert first_focusing_quad_mixed.k == pytest.approx(
        1.32
    )  # Modified value from type-based filtering


def test_performance_with_large_lattice():
    """Test performance and correctness with a larger lattice."""

    # Create a larger lattice
    lattice = elements.KnownElementsList()

    # Add many elements
    for i in range(100):
        lattice.append(elements.Drift(name=f"drift{i}", ds=1.0))
        lattice.append(
            elements.Quad(name=f"quad{i}", ds=0.3, k=1.0 if i % 2 == 0 else -1.0)
        )
        lattice.append(elements.Marker(name=f"marker{i}"))

    # Test filtering performance
    drift_elements = lattice.select(kind="Drift")
    assert len(drift_elements) == 100

    quad_elements = lattice.select(kind="Quad")
    assert len(quad_elements) == 100

    marker_elements = lattice.select(kind="Marker")
    assert len(marker_elements) == 100

    # Test OR filtering
    drift_quad_elements = lattice.select(kind=["Drift", "Quad"])
    assert len(drift_quad_elements) == 200

    # Test reference preservation in large lattice
    first_drift = lattice.select(kind="Drift")[0]
    first_drift.ds = 999.0
    assert lattice[0].ds == 999.0  # Original element modified


def test_regex_name_filtering():
    """Test regex pattern filtering for element names."""

    # Create a test lattice with various naming patterns
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Quad(name="quad2", ds=0.3, k=-1.0),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
            elements.Quad(name="qf1", ds=0.2, k=1.5),
            elements.Quad(name="qd1", ds=0.2, k=-1.5),
            elements.Marker(name="marker1"),
            elements.Drift(name="drift_long_name", ds=3.0),
        ]
    )

    # Test single regex pattern
    quad_elements = lattice.select(name=r"quad\d+")
    assert len(quad_elements) == 2
    assert all(el.name in ["quad1", "quad2"] for el in quad_elements)

    # Test regex pattern for drift elements
    drift_pattern = r"drift\d+"
    drift_elements = lattice.select(name=drift_pattern)
    assert len(drift_elements) == 2
    assert all(el.name in ["drift1", "drift2"] for el in drift_elements)

    # Test regex pattern that matches multiple patterns
    q_pattern = r"q[fd]\d+"
    q_elements = lattice.select(name=q_pattern)
    assert len(q_elements) == 2
    assert all(el.name in ["qf1", "qd1"] for el in q_elements)

    # Test regex pattern with no matches
    no_match_pattern = r"nonexistent\d+"
    no_match_elements = lattice.select(name=no_match_pattern)
    assert len(no_match_elements) == 0

    # Test regex pattern with word boundaries
    exact_pattern = r"\bdrift\d+\b"
    exact_drift_elements = lattice.select(name=exact_pattern)
    assert len(exact_drift_elements) == 2
    assert all(el.name in ["drift1", "drift2"] for el in exact_drift_elements)

    # Test regex pattern that should not match "drift_long_name"
    short_pattern = r"drift\d$"
    short_drift_elements = lattice.select(name=short_pattern)
    assert len(short_drift_elements) == 2
    assert all(el.name in ["drift1", "drift2"] for el in short_drift_elements)


def test_mixed_regex_string_filtering():
    """Test mixing regex patterns and strings in name filtering."""

    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.append(elements.Drift(name="drift1", ds=1.0))
    lattice.append(elements.Quad(name="quad1", ds=0.5, k=1.0))
    lattice.append(elements.Drift(name="drift2", ds=2.0))
    lattice.append(elements.Quad(name="quad2", ds=0.3, k=-1.0))
    lattice.append(elements.Sbend(name="bend1", ds=1.0, rc=10.0))
    lattice.append(elements.Quad(name="qf1", ds=0.2, k=1.5))
    lattice.append(elements.Marker(name="marker1"))

    # Test mixing regex and string patterns
    mixed_patterns = [r"quad\d+", "bend1"]
    mixed_elements = lattice.select(name=mixed_patterns)
    assert len(mixed_elements) == 3
    assert all(el.name in ["quad1", "quad2", "bend1"] for el in mixed_elements)

    # Test mixing regex and string patterns with tuple
    mixed_tuple = (r"drift\d+", "marker1")
    mixed_tuple_elements = lattice.select(name=mixed_tuple)
    assert len(mixed_tuple_elements) == 3
    assert all(
        el.name in ["drift1", "drift2", "marker1"] for el in mixed_tuple_elements
    )

    # Test multiple regex patterns
    multi_regex = [r"quad\d+", r"qf\d+"]
    multi_regex_elements = lattice.select(name=multi_regex)
    assert len(multi_regex_elements) == 3
    assert all(el.name in ["quad1", "quad2", "qf1"] for el in multi_regex_elements)


def test_regex_error_handling():
    """Test error handling for invalid regex usage."""

    lattice = elements.KnownElementsList()
    lattice.append(elements.Drift(name="drift1", ds=1.0))

    # Test invalid name parameter type
    with pytest.raises(
        TypeError,
        match="'name' parameter must be a string or list/tuple of strings",
    ):
        _ = lattice.select(name=123)

    # Test invalid name parameter with mixed types in list
    with pytest.raises(TypeError, match="'name' parameter must contain strings"):
        _ = lattice.select(name=["drift1", 123])

    # Test invalid name parameter with mixed types in list including regex
    with pytest.raises(TypeError, match="'name' parameter must contain strings"):
        _ = lattice.select(name=[r"drift\d+", 123])


def test_regex_chaining():
    """Test chaining filters with regex patterns."""

    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Quad(name="quad2", ds=0.3, k=-1.0),
            elements.Sbend(name="bend1", ds=1.0, rc=10.0),
            elements.Quad(name="qf1", ds=0.2, k=1.5),
        ]
    )

    # Test chaining kind filter then regex name filter
    quad_then_regex = lattice.select(kind="Quad").select(name=r"quad\d+")
    assert len(quad_then_regex) == 2
    assert all(el.name in ["quad1", "quad2"] for el in quad_then_regex)
    assert all(type(el) is elements.Quad for el in quad_then_regex)

    # Test chaining regex name filter then kind filter
    regex_then_kind = lattice.select(name=r"q\w+\d+").select(kind="Quad")
    assert len(regex_then_kind) == 3
    assert all(el.name in ["quad1", "quad2", "qf1"] for el in regex_then_kind)
    assert all(type(el) is elements.Quad for el in regex_then_kind)

    # Test chaining with mixed patterns
    mixed_then_kind = lattice.select(name=[r"quad\d+", "bend1"]).select(kind="Quad")
    assert len(mixed_then_kind) == 2
    assert all(el.name in ["quad1", "quad2"] for el in mixed_then_kind)
    assert all(type(el) is elements.Quad for el in mixed_then_kind)

    # Test reference preservation for regex filtering
    quad_regex = lattice.select(name=r"quad\d+")
    quad_regex[0].k = 2.0
    assert lattice[1].k == 2.0  # Original element should be modified

    # Test reference preservation for mixed patterns
    mixed_patterns = lattice.select(name=[r"drift\d+", "quad2"])
    mixed_patterns[0].ds = 5.0
    assert lattice[0].ds == 5.0  # Original element should be modified

    # Test that modifying through different regex patterns affects the same element
    drift_regex = lattice.select(name=r"drift\d+")
    drift_regex[0].ds = 10.0
    assert lattice[0].ds == 10.0  # Same element modified through different filter

    # Test manipulation of chained regex results
    quad_then_regex[0].k = 5.0  # Modify quad1
    quad_then_regex[1].k = -5.0  # Modify quad2
    assert lattice[1].k == 5.0  # quad1 should be modified
    assert lattice[3].k == -5.0  # quad2 should be modified

    # Test manipulation of regex_then_kind results
    regex_then_kind[2].k = 3.0  # Modify qf1
    assert lattice[5].k == 3.0  # qf1 should be modified


def test_select_no_arguments():
    """Test select() method with no arguments returns all elements."""
    # Create a test lattice
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Quad(name="quad2", ds=0.3, k=-1.0),
        ]
    )

    # Test select() with no arguments returns all elements
    all_elements = lattice.select()
    assert len(all_elements) == 4
    assert [el.name for el in all_elements] == ["drift1", "quad1", "drift2", "quad2"]

    # Test modifications affect original elements
    all_elements[0].ds = 1.5
    assert lattice[0].ds == 1.5
    lattice[0].ds = 1.0  # Reset

    # Test chaining: select(something).select() and select().select(something)
    drift_then_all = lattice.select(kind="Drift").select()
    assert len(drift_then_all) == 2  # Still only drifts
    assert [el.name for el in drift_then_all] == ["drift1", "drift2"]

    all_then_drift = lattice.select().select(kind="Drift")
    assert len(all_then_drift) == 2
    assert [el.name for el in all_then_drift] == ["drift1", "drift2"]


def test_delete_removes_and_returns_none():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
        ]
    )
    sel = lattice.select(kind="Drift")
    ret = sel.delete()
    assert ret is None
    assert len(lattice) == 1
    assert type(lattice[0]) is elements.Quad
    assert lattice[0].name == "quad1"


def test_delete_invalidates_chained_filtered_views():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
        ]
    )
    outer = lattice.select(kind="Drift")
    inner = outer.select(name="drift1")
    inner.delete()
    with pytest.raises(RuntimeError, match="no longer valid"):
        _ = len(outer)
    with pytest.raises(RuntimeError, match="no longer valid"):
        _ = len(inner)
    assert len(lattice) == 2
    assert lattice[0].name == "quad1"
    assert lattice[1].name == "drift2"


def test_delete_invalidates_sibling_filtered_views():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
        ]
    )
    drifts = lattice.select(kind="Drift")
    branch_a = drifts.select(name="drift1")
    branch_b = drifts.select(name="drift2")
    branch_a.delete()
    with pytest.raises(RuntimeError, match="no longer valid"):
        _ = len(drifts)
    with pytest.raises(RuntimeError, match="no longer valid"):
        _ = len(branch_b)
    assert len(lattice) == 2


def test_replace_each_default_keep_name_not_ds():
    """Default: keep_name=True, keep_ds=False

    names from replaced elements, ds from template."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="qf1", ds=0.5, k=1.0),
            elements.Quad(name="qd1", ds=0.25, k=-1.0),
        ]
    )
    template = elements.Quad(name="tpl", ds=0.01, k=99.0)
    sel = lattice.select(kind="Quad")
    new_sel = sel.replace_each(template)
    assert type(lattice[1]) is elements.Quad
    assert lattice[1].name == "qf1"
    assert lattice[1].ds == pytest.approx(0.01)
    assert lattice[1].k == 99.0
    assert lattice[2].name == "qd1"
    assert lattice[2].ds == pytest.approx(0.01)
    assert len(new_sel) == 2
    assert new_sel[0].name == "qf1"
    with pytest.raises(RuntimeError, match="no longer valid"):
        _ = len(sel)


def test_replace_each_keep_ds_true():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(name="qf1", ds=0.5, k=1.0),
            elements.Quad(name="qd1", ds=0.25, k=-1.0),
        ]
    )
    template = elements.Quad(name="tpl", ds=0.01, k=99.0)
    lattice.select(kind="Quad").replace_each(template, keep_ds=True)
    assert lattice[0].ds == 0.5
    assert lattice[1].ds == 0.25
    assert lattice[0].k == 99.0


def test_replace_each_keep_name_and_ds_false():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(name="qf1", ds=0.5, k=1.0),
        ]
    )
    template = elements.Quad(name="tpl", ds=0.01, k=7.0)
    lattice.select(kind="Quad").replace_each(template, keep_name=False)
    assert lattice[0].name == "tpl"
    assert lattice[0].ds == pytest.approx(0.01)
    assert lattice[0].k == 7.0


def test_replace_with_drifts_model_match():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(name="lq", ds=1.0, k=1.0, nslice=2),
            elements.ChrQuad(name="cq", ds=0.5, k=1.0, nslice=3),
            elements.ExactQuad(name="eq", ds=0.25, k=1.0, nslice=4),
        ]
    )
    sel = lattice.select(kind=r".*Quad")
    sel.replace_with_drifts()
    assert type(lattice[0]) is elements.Drift
    assert type(lattice[1]) is elements.ChrDrift
    assert type(lattice[2]) is elements.ExactDrift
    assert lattice[0].name == "lq"
    assert lattice[1].name == "cq"
    assert lattice[2].name == "eq"
    assert lattice[0].ds == 1.0
    assert lattice[1].ds == 0.5
    assert lattice[2].ds == 0.25
    assert lattice[0].nslice == 2
    assert lattice[1].nslice == 3
    assert lattice[2].nslice == 4


def test_replace_with_drifts_explicit_model():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(name="lq", ds=1.0, k=1.0),
            elements.ExactQuad(name="eq", ds=0.25, k=1.0),
        ]
    )
    lattice.select(kind=r".*Quad").replace_with_drifts(model="paraxial")
    assert type(lattice[0]) is elements.ChrDrift
    assert type(lattice[1]) is elements.ChrDrift
    assert lattice[0].ds == 1.0
    assert lattice[1].ds == 0.25

    lattice.select(kind=r".*Drift").replace_with_drifts(model="exact")
    assert type(lattice[0]) is elements.ExactDrift
    assert type(lattice[1]) is elements.ExactDrift


def test_replace_invalid_model_string():
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Drift(name="d1", ds=1.0)])
    sel = lattice.select(kind="Drift")
    with pytest.raises(ValueError, match="model"):
        sel.replace_with_drifts(model="not_a_model")


def test_replace_with_drifts_keep_alignment():
    """Test keep_alignment parameter."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [elements.Quad(name="q1", ds=0.5, k=1.0, dx=0.01, dy=0.02, rotation=0.5)]
    )

    # Default: keep_alignment=True
    lattice.select(kind="Quad").replace_with_drifts()
    assert lattice[0].dx == pytest.approx(0.01)
    assert lattice[0].dy == pytest.approx(0.02)
    assert lattice[0].rotation == 0.5

    # Reset and test keep_alignment=False
    lattice.clear()
    lattice.extend(
        [elements.Quad(name="q1", ds=0.5, k=1.0, dx=0.01, dy=0.02, rotation=0.5)]
    )
    lattice.select(kind="Quad").replace_with_drifts(keep_alignment=False)
    assert lattice[0].dx == 0.0
    assert lattice[0].dy == 0.0
    assert lattice[0].rotation == 0.0


def test_replace_with_drifts_keep_aperture():
    """Test keep_aperture parameter."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [elements.Quad(name="q1", ds=0.5, k=1.0, aperture_x=0.05, aperture_y=0.03)]
    )

    # Default: keep_aperture=False
    lattice.select(kind="Quad").replace_with_drifts()
    assert lattice[0].aperture_x == 0.0
    assert lattice[0].aperture_y == 0.0

    # Reset and test keep_aperture=True
    lattice.clear()
    lattice.extend(
        [elements.Quad(name="q1", ds=0.5, k=1.0, aperture_x=0.05, aperture_y=0.03)]
    )
    lattice.select(kind="Quad").replace_with_drifts(keep_aperture=True)
    assert lattice[0].aperture_x == pytest.approx(0.05)
    assert lattice[0].aperture_y == pytest.approx(0.03)


def test_chained_replace_each():
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(name="qf1", ds=0.5, k=1.0),
            elements.Quad(name="qd1", ds=0.25, k=-1.0),
        ]
    )
    q2 = elements.Quad(name="tpl", ds=0.1, k=2.0)
    new_sel = (
        lattice.select(kind="Quad")
        .replace_each(q2)
        .replace_each(elements.Quad(name="tpl2", ds=0.2, k=3.0))
    )
    assert lattice[0].k == 3.0
    assert lattice[1].k == 3.0
    assert len(new_sel) == 2


def test_delete_empty_selection():
    """delete() on empty selection is a no-op, returns None."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="d1", ds=1.0),
            elements.Quad(name="q1", ds=0.5, k=1.0),
        ]
    )
    sel = lattice.select(name="nonexistent")
    assert len(sel) == 0
    ret = sel.delete()
    assert ret is None
    assert len(lattice) == 2


def test_replace_each_empty_selection():
    """replace_each on empty selection returns empty FilteredElementsList."""
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Drift(name="d1", ds=1.0)])
    sel = lattice.select(name="nonexistent")
    new_sel = sel.replace_each(elements.Quad(name="tpl", ds=0.1, k=1.0))
    assert len(new_sel) == 0
    assert len(lattice) == 1
    assert lattice[0].name == "d1"


def test_replace_with_drifts_empty_selection():
    """replace_with_drifts on empty selection returns empty FilteredElementsList."""
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Drift(name="d1", ds=1.0)])
    sel = lattice.select(name="nonexistent")
    new_sel = sel.replace_with_drifts()
    assert len(new_sel) == 0
    assert len(lattice) == 1


def test_replace_each_thin_element():
    """replace_each with keep_ds=True on thin element (Marker has ds=0)."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Marker(name="m1"),
            elements.Drift(name="d1", ds=1.0),
        ]
    )
    # Marker is a thin element with ds=0; keep_ds=True should preserve that
    sel = lattice.select(kind="Marker")
    new_sel = sel.replace_each(elements.Drift(name="tpl", ds=0.5), keep_ds=True)
    assert len(new_sel) == 1
    assert type(lattice[0]) is elements.Drift
    assert lattice[0].name == "m1"
    # Marker has ds=0.0, so with keep_ds=True, the new element gets ds=0.0
    assert lattice[0].ds == 0.0

    # Now test keep_ds=False: template's ds should be used
    lattice2 = elements.KnownElementsList()
    lattice2.extend([elements.Marker(name="m2")])
    lattice2.select(kind="Marker").replace_each(
        elements.Drift(name="tpl", ds=0.5), keep_ds=False
    )
    assert lattice2[0].ds == 0.5


def test_replace_with_drifts_thin_element():
    """replace_with_drifts on thin element (Marker with ds=0)."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Marker(name="m1"),
            elements.Drift(name="d1", ds=1.0),
        ]
    )
    sel = lattice.select(kind="Marker")
    new_sel = sel.replace_with_drifts()
    assert len(new_sel) == 1
    assert type(lattice[0]) is elements.Drift
    assert lattice[0].name == "m1"
    # Marker has ds=0.0 (thin element), so drift gets ds=0.0
    assert lattice[0].ds == 0.0
    assert lattice[1].name == "d1"
    assert lattice[1].ds == 1.0


def test_delete_with_thin_element_present():
    """Cloning must not pass `ds` to thin elements that do not accept it.

    `Marker.to_dict()` reports `ds=0.0` but `Marker.__init__` takes only a name,
    so rebuilding the lattice used to raise TypeError whenever an unselected
    thin element had to be cloned.
    """
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Marker(name="m1"),
            elements.Drift(name="d1", ds=1.0),
            elements.ShortRF(name="rf1", V=1.0, freq=1.0e6, phase=0.0),
        ]
    )
    lattice.select(kind="Drift").delete()
    assert [type(e).__name__ for e in lattice] == ["Marker", "ShortRF"]
    assert lattice[0].name == "m1"
    assert lattice[1].name == "rf1"


def test_replace_each_with_thin_element_present():
    """Same cloning path, via replace_each."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Marker(name="m1"),
            elements.Quad(name="q1", ds=0.5, k=1.0),
        ]
    )
    lattice.select(kind="Quad").replace_each(elements.Drift(ds=0.5))
    assert [type(e).__name__ for e in lattice] == ["Marker", "Drift"]
    assert lattice[0].name == "m1"


# ---------------------------------------------------------------------------
# select(has=...) capability filter
# ---------------------------------------------------------------------------


def _mixed_lattice():
    """Thick, exact and thin elements, covering every settability case."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="d1", ds=1.0),
            elements.Quad(name="q1", ds=0.5, k=1.0),
            elements.ExactQuad(name="eq1", ds=0.5, k=1.0),
            elements.ExactDrift(name="ed1", ds=0.3),
            elements.Marker(name="m1"),
            elements.ShortRF(name="rf1", V=1.0, freq=1.0e6, phase=0.0),
        ]
    )
    return lattice


def _kinds(selection):
    return [type(e).__name__ for e in selection]


def test_select_has_is_settable_not_merely_readable():
    """`nslice` reads on every element but only thick ones accept a value."""
    lattice = _mixed_lattice()
    # every element reports nslice ...
    assert all(hasattr(e, "nslice") for e in lattice)
    # ... but has= matches only those that can be assigned
    assert _kinds(lattice.select(has="nslice")) == [
        "Drift",
        "Quad",
        "ExactQuad",
        "ExactDrift",
    ]


def test_select_has_list_is_or_matched():
    lattice = _mixed_lattice()
    lattice.append(
        elements.SoftQuadrupole(
            name="sq1",
            ds=1.0,
            gscale=1.0,
            cos_coefficients=[1.0],
            sin_coefficients=[0.0],
        )
    )
    # SoftQuadrupole has mapsteps but no int_order; ExactQuad has both
    assert _kinds(lattice.select(has="mapsteps")) == ["ExactQuad", "SoftQuadrupole"]
    assert _kinds(lattice.select(has=["int_order", "mapsteps"])) == [
        "ExactQuad",
        "SoftQuadrupole",
    ]


def test_select_has_ors_with_kind_and_ands_when_chained():
    lattice = _mixed_lattice()
    # OR within a single call, like kind and name
    assert _kinds(lattice.select(kind="Marker", has="int_order")) == [
        "ExactQuad",
        "Marker",
    ]
    # AND by chaining
    assert _kinds(lattice.select(kind=r".*Quad").select(has="int_order")) == [
        "ExactQuad"
    ]


def test_select_has_invalid_type_raises():
    lattice = _mixed_lattice()
    with pytest.raises(TypeError, match="has"):
        lattice.select(has=42)
    with pytest.raises(TypeError, match="has"):
        lattice.select(has=["nslice", 42])


# ---------------------------------------------------------------------------
# set()
# ---------------------------------------------------------------------------


def test_set_on_selection():
    lattice = _mixed_lattice()
    assert lattice.select(kind=r".*Quad").set(nslice=8) == 2
    assert [e.nslice for e in lattice] == [1, 8, 8, 1, 1, 1]


def test_set_raises_on_mixed_lattice_and_changes_nothing():
    """Default is raise, not skip: a whole-lattice call must not half-apply."""
    lattice = _mixed_lattice()
    with pytest.raises(AttributeError, match="nslice"):
        lattice.set(nslice=8)
    assert [e.nslice for e in lattice] == [1, 1, 1, 1, 1, 1]


def test_set_skip_applies_where_possible():
    lattice = _mixed_lattice()
    assert lattice.set(nslice=8, skip=True) == 4
    assert [e.nslice for e in lattice] == [8, 8, 8, 8, 1, 1]


def test_set_multiple_properties():
    lattice = _mixed_lattice()
    assert lattice.select(has="int_order").set(int_order=6, mapsteps=6) == 1
    eq = lattice[2]
    assert eq.int_order == 6
    assert eq.mapsteps == 6


def test_set_unknown_property_raises_by_default():
    """A typo names a property no element has, which must not pass silently."""
    lattice = _mixed_lattice()
    with pytest.raises(AttributeError, match="nslize"):
        lattice.set(nslize=8)
    assert lattice.set(nslize=8, skip=True) == 0


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (dict(nslice=0), "nslice"),
        (dict(nslice=-3), "nslice"),
        (dict(int_order=3), "int_order"),
        (dict(mapsteps=0), "mapsteps"),
    ],
)
def test_set_rejects_invalid_values_without_writing(kwargs, match):
    lattice = _mixed_lattice()
    before = [(e.nslice, getattr(e, "int_order", None)) for e in lattice]
    with pytest.raises(ValueError, match=match):
        lattice.select(has="int_order").set(**kwargs)
    assert [(e.nslice, getattr(e, "int_order", None)) for e in lattice] == before


def test_set_does_not_invalidate_live_selections():
    """Contrast with delete(), which rebuilds the lattice and invalidates views."""
    lattice = _mixed_lattice()
    view = lattice.select(kind="Quad")
    lattice.set(nslice=4, skip=True)
    # the view is still usable and sees the new value
    assert len(view) == 1
    assert view[0].nslice == 4

    lattice.select(kind="Drift").delete()
    with pytest.raises(RuntimeError, match="no longer valid"):
        len(view)


def test_set_on_empty_selection():
    lattice = _mixed_lattice()
    empty = lattice.select(name="does_not_exist")
    assert len(empty) == 0
    assert empty.set(nslice=8) == 0


def test_set_via_filtered_view_matches_container():
    """set() is registered on both containers with the same semantics."""
    lattice = _mixed_lattice()
    assert lattice.select().set(nslice=2, skip=True) == 4
    assert lattice.set(nslice=3, skip=True) == 4
    assert [e.nslice for e in lattice] == [3, 3, 3, 3, 1, 1]


# ---------------------------------------------------------------------------
# element-level validation these rely on
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "ctor",
    [
        lambda: elements.Drift(ds=1.0, nslice=0),
        lambda: elements.ExactQuad(ds=1.0, k=1.0, int_order=3),
        lambda: elements.ExactQuad(ds=1.0, k=1.0, mapsteps=0),
        lambda: elements.ExactMultipole(
            ds=1.0, k_normal=[0.1], k_skew=[0.0], int_order=5
        ),
    ],
)
def test_constructors_reject_invalid_knobs(ctor):
    with pytest.raises(ValueError):
        ctor()


def test_property_setters_reject_invalid_knobs():
    d = elements.Drift(ds=1.0)
    with pytest.raises(ValueError):
        d.nslice = 0
    assert d.nslice == 1

    q = elements.ExactQuad(ds=1.0, k=1.0)
    with pytest.raises(ValueError):
        q.int_order = 3
    with pytest.raises(ValueError):
        q.mapsteps = 0
    assert (q.int_order, q.mapsteps) == (2, 10)


def test_argument_validation_raises_value_error():
    """Element argument checks raise ValueError, not RuntimeError."""
    with pytest.raises(ValueError, match="R must be > 0"):
        elements.DipEdge(psi=0.1, rc=1.0, g=0.01, R=-1.0)
    with pytest.raises(ValueError, match="same length"):
        elements.ExactCFbend(ds=1.0, k_normal=[0.1], k_skew=[0.0, 0.0])


def test_set_rolls_back_when_a_setter_rejects_late():
    """A setter that rejects after earlier writes must undo them.

    `_SET_VALIDATORS` only prevalidates nslice/int_order/mapsteps. `DipEdge.R`
    is validated by its own setter, so `set(psi=..., R=-1.0)` fails only after
    `psi` has been written. That write has to be rolled back.
    """
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.DipEdge(psi=0.1, rc=1.0, g=0.01, R=1.0, name="de1"),
            elements.DipEdge(psi=0.1, rc=1.0, g=0.01, R=1.0, name="de2"),
        ]
    )
    before = [(e.psi, e.R) for e in lattice]

    with pytest.raises(ValueError, match="R must be > 0"):
        lattice.set(psi=0.2, R=-1.0)

    assert [(e.psi, e.R) for e in lattice] == before


def test_set_rolls_back_earlier_elements_too():
    """The rollback spans elements, not just properties within one element."""
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            # a first element that accepts everything
            elements.DipEdge(psi=0.1, rc=1.0, g=0.01, R=1.0, name="ok"),
            # a second one that will reject R and must undo the first
            elements.DipEdge(psi=0.1, rc=1.0, g=0.01, R=1.0, name="bad"),
        ]
    )
    with pytest.raises(ValueError):
        lattice.set(psi=0.3, R=-2.0)

    assert [e.psi for e in lattice] == [0.1, 0.1]


def test_lattice_stores_copies_not_references():
    """Documented container semantics that `set()` inherits.

    `append` copies, so the lattice element and the object handed to `append`
    are independent. A `Programmable` callback closing over the original will
    therefore not see `ds`/`nslice` written through the lattice. Tracked in
    https://github.com/BLAST-ImpactX/impactx/pull/1381
    """
    original = elements.Drift(ds=1.0, name="d")
    lattice = elements.KnownElementsList()
    lattice.append(original)

    assert lattice[0] is not original
    assert lattice[0] is lattice[0]  # but stable within the lattice

    lattice.set(ds=7.0)
    assert lattice[0].ds == 7.0
    assert original.ds == 1.0  # unchanged

    original.ds = 9.0
    assert lattice[0].ds == 7.0  # and not the reverse either
