#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""Array-valued element parameters: Fourier and multipole coefficients, polygon vertices.

They are set after construction through a property each, or together through a paired
setter, validated before anything is stored, and reach the next particle push.
"""

import numpy as np
import pytest

from impactx import ImpactX, RefPart, elements

# element name -> (settings besides the arrays, first array, second array)
COEFFICIENT_ELEMENTS = {
    "SoftQuadrupole": ({"gscale": 1.0}, "cos_coefficients", "sin_coefficients"),
    "SoftSolenoid": ({"bscale": 1.0}, "cos_coefficients", "sin_coefficients"),
    "RFCavity": (
        {"escale": 1.0, "freq": 1.0e9, "phase": 10.0},
        "cos_coefficients",
        "sin_coefficients",
    ),
    "ExactMultipole": ({}, "k_normal", "k_skew"),
    "ExactCFbend": ({}, "k_normal", "k_skew"),
}
NAMES = sorted(COEFFICIENT_ELEMENTS)


def build(name, first, second):
    settings, first_key, second_key = COEFFICIENT_ELEMENTS[name]
    return getattr(elements, name)(
        ds=0.1, **settings, **{first_key: first, second_key: second}
    )


@pytest.mark.parametrize("name", NAMES)
def test_coefficients_are_set_and_validated(name):
    _, first, second = COEFFICIENT_ELEMENTS[name]
    el = build(name, [1.0, 2.0], [0.0, 3.0])
    assert (getattr(el, first), getattr(el, second)) == ([1.0, 2.0], [0.0, 3.0])

    # one array at a time keeps the other
    setattr(el, first, [5.0, 6.0])
    assert (getattr(el, first), getattr(el, second)) == ([5.0, 6.0], [0.0, 3.0])

    # the paired setter can change the length
    el.set_coefficients([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
    assert getattr(el, first) == [1.0, 2.0, 3.0]

    # rejected updates change nothing
    with pytest.raises(ValueError, match="same length"):
        setattr(el, first, [1.0])
    with pytest.raises(ValueError, match="same length"):
        el.set_coefficients([1.0, 2.0], [1.0])
    if name != "ExactMultipole":  # the field evaluation reads the first coefficient
        with pytest.raises(ValueError):
            el.set_coefficients([], [])
        with pytest.raises(ValueError):
            build(name, [], [])
    assert getattr(el, first) == [1.0, 2.0, 3.0]
    assert getattr(el, second) == [0.0, 0.0, 0.0]


def test_a_single_multipole_coefficient_has_a_transfer_map():
    """A pure dipole supplies no quadrupole coefficient."""

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    elements.ExactMultipole(ds=0.1, k_normal=[1.0], k_skew=[0.0]).transfer_map(ref)


def test_polygon_vertices_are_set_and_validated():
    outline_x = [-1.0, 1.0, 1.0, -1.0, -1.0]
    outline_y = [-1.0, -1.0, 1.0, 1.0, -1.0]
    poly = elements.PolygonAperture(vertices_x=outline_x, vertices_y=outline_y)
    assert (poly.vertices_x, poly.vertices_y) == (outline_x, outline_y)

    with pytest.raises(ValueError, match="same length"):
        poly.set_vertices(outline_x, [])
    with pytest.raises(ValueError, match="same length"):
        poly.vertices_y = [0.0, 1.0]
    with pytest.raises(ValueError, match="first and last vertex"):
        poly.set_vertices([0.0, 1.0, 2.0], [0.0, 1.0, 2.0])
    assert (poly.vertices_x, poly.vertices_y) == (outline_x, outline_y)

    poly.set_vertices([0.0, 2.0, 2.0, 0.0, 0.0], [0.0, 0.0, 2.0, 2.0, 0.0])
    assert poly.vertices_x == [0.0, 2.0, 2.0, 0.0, 0.0]


@pytest.fixture
def push_once():
    """Push three fixed particles through an element and return the result."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.n_cell = [8, 8, 8]
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    def run(element):
        beam = sim.beam
        beam.clear_particles()
        beam.ref.reset()
        beam.ref.set_species("electron").set_kin_energy_MeV(100.0)
        beam.add_n_particles(
            [0.0, 1.0e-3, 3.0e-3],
            [0.0, 2.0e-3, -1.0e-3],
            [0.0, 3.0e-3, -2.0e-3],
            [1.0e-4, 2.0e-4, -1.0e-4],
            [2.0e-4, -1.0e-4, 3.0e-4],
            [3.0e-4, 1.0e-4, -2.0e-4],
            -1.0 / 0.510998950e6,
            1.0e-12,
        )
        element.push(beam)
        phase_space = beam.to_df(local=True)[
            [
                "position_x",
                "position_y",
                "position_t",
                "momentum_x",
                "momentum_y",
                "momentum_t",
            ]
        ].to_numpy(copy=True)
        alive = beam.total_number_of_particles(only_valid=True, only_local=True)
        return phase_space, alive

    try:
        yield run
    finally:
        sim.finalize()


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize(
    "first,second",
    [([2.0, 3.0], [0.0, 0.1]), ([2.0, 3.0, 0.4], [0.0, 0.1, 0.2])],
    ids=["same_length", "grow"],
)
def test_new_coefficients_reach_the_particle_push(name, first, second, push_once):
    """An element pushed before and after an update pushes like a fresh element."""

    element = build(name, [1.0, 2.0], [0.0, 3.0])
    before, _ = push_once(element)

    element.set_coefficients(first, second)
    after, alive = push_once(element)
    expected, _ = push_once(build(name, first, second))

    assert alive == 3
    assert np.isfinite(after).all()
    assert not np.array_equal(before, after)
    # identical arithmetic on the same backend
    np.testing.assert_array_equal(after, expected)


def test_new_vertices_reach_the_particle_push(push_once):
    polygon = elements.PolygonAperture(
        vertices_x=[-0.01, 0.01, 0.01, -0.01, -0.01],
        vertices_y=[-0.01, -0.01, 0.01, 0.01, -0.01],
    )
    _, before = push_once(polygon)
    assert before == 3

    narrow_x = [-0.0015, 0.0, 0.0015, 0.0015, -0.0015, -0.0015]
    narrow_y = [-0.004, -0.004, -0.004, 0.004, 0.004, -0.004]
    polygon.set_vertices(narrow_x, narrow_y)
    _, after = push_once(polygon)

    assert after == 2
