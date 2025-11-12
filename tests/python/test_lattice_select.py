#!/usr/bin/env python3
#
# Copyright 2022-2024 The ImpactX Community
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
    assert lattice[1].k == 1.1  # qf1
    assert lattice[5].k == -1.1  # qd1
    assert lattice[7].k == 1.1  # qf2

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
    assert lattice[3].ds == 0.9  # bend1
    assert lattice[9].ds == 0.9  # bend2

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


def test_position_filtering():
    """Test position-based filtering with s parameter."""
    import impactx
    from impactx import elements

    # Create a lattice with known positions
    lattice = impactx.elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="drift1", ds=1.0),
            elements.Quad(name="quad1", ds=0.5, k=1.0),
            elements.Drift(name="drift2", ds=2.0),
            elements.Quad(name="quad2", ds=0.3, k=-1.0),
            elements.Drift(name="drift3", ds=1.5),
        ]
    )

    # Test range filtering
    early_elements = lattice.select(s=(None, 2.0))
    assert len(early_elements) == 3  # drift1, quad1, drift2
    assert early_elements[0].name == "drift1"
    assert early_elements[1].name == "quad1"
    assert early_elements[2].name == "drift2"

    # Test range filtering with both bounds (overlap logic)
    middle_elements = lattice.select(s=(1.0, 3.0))
    assert len(middle_elements) == 3  # drift1, quad1, drift2 (all overlap with range)
    assert middle_elements[0].name == "drift1"  # s=[0.0, 1.0] overlaps with [1.0, 3.0]
    assert middle_elements[1].name == "quad1"  # s=[1.0, 1.5] overlaps with [1.0, 3.0]
    assert middle_elements[2].name == "drift2"  # s=[1.5, 3.5] overlaps with [1.0, 3.0]

    # Test upper bound only
    late_elements = lattice.select(s=(3.0, None))
    assert len(late_elements) == 3  # drift2, quad2, drift3
    assert late_elements[0].name == "drift2"
    assert late_elements[1].name == "quad2"
    assert late_elements[2].name == "drift3"

    # Test combined filtering (OR logic)
    drift_or_early = lattice.select(kind="Drift", s=(None, 2.0))
    assert len(drift_or_early) == 4  # All drifts + early elements
    drift_names = [el.name for el in drift_or_early if el.name.startswith("drift")]
    assert len(drift_names) == 3  # All 3 drifts
    assert "drift1" in drift_names
    assert "drift2" in drift_names
    assert "drift3" in drift_names

    # Test chaining (AND logic)
    drift_then_early = lattice.select(kind="Drift").select(s=(None, 1.9))
    assert len(drift_then_early) == 2  # Only drift1 and drift2
    assert drift_then_early[0].name == "drift1"
    assert drift_then_early[1].name == "drift2"

    # Test reference preservation through chaining
    original_ds = lattice[1].ds
    chained_elements = lattice.select(kind="Quad").select(s=(1.0, 1.9))
    chained_elements[0].ds = 0.8  # Modify through chained filter
    assert lattice[1].ds == 0.8  # Original element modified
    lattice[1].ds = original_ds  # Reset

    # Test with list syntax
    range_list = lattice.select(s=[1.0, 2.9])
    assert len(range_list) == 3
    assert range_list[0].name == "drift1"

    # Test explicit overlap behavior
    # Element at s=[0.0, 1.0] should be included in range [0.5, 1.5] because it overlaps
    # Element at s=[1.0, 1.5] should be included in range [0.5, 1.5] because it overlaps
    # Element at s=[1.5, 3.5] should NOT be included because it doesn't overlap with [0.5, 1.5]
    overlap_test = lattice.select(s=(0.5, 1.4))
    assert len(overlap_test) == 2  # drift1 and quad1 both overlap
    assert overlap_test[0].name == "drift1"  # s=[0.0, 1.0] overlaps with [0.5, 1.5]
    assert overlap_test[1].name == "quad1"  # s=[1.0, 1.5] overlaps with [0.5, 1.5]

    # Test non-overlap case
    # Range [0.2, 0.8] should only include drift1 (s=[0.0, 1.0])
    narrow_test = lattice.select(s=(0.2, 0.8))
    assert len(narrow_test) == 1
    assert narrow_test[0].name == "drift1"  # Only drift1 overlaps with [0.2, 0.8]

    # Test error handling
    with pytest.raises(TypeError, match="must have exactly 2 elements"):
        lattice.select(s=(1.0,))  # Wrong number of elements

    with pytest.raises(TypeError, match="must be a tuple/list with 2 elements"):
        lattice.select(s="invalid")  # Wrong type

    with pytest.raises(TypeError, match="bounds must be numbers or None"):
        lattice.select(s=(1.0, "invalid"))  # Wrong bound type
