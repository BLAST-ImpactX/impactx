#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""Survey of a lattice: reference trajectory and envelope."""

import math

import numpy as np
import pandas as pd
import pytest

from impactx import Config, ImpactX, distribution, elements, push

SIGMAS = ("sigma_x", "sigma_px", "sigma_y", "sigma_py", "sigma_t", "sigma_pt")


def make_sim(npart=10000):
    """A 100 MeV electron beam with correlations in each plane."""
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.verbose = 0
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(100.0)
    distr = distribution.Waterbag(
        lambdaX=1.0e-4,
        lambdaY=2.0e-4,
        lambdaT=1.0e-3,
        lambdaPx=3.0e-5,
        lambdaPy=1.0e-5,
        lambdaPt=2.0e-3,
        muxpx=-0.5,
        muypy=0.3,
        mutpt=0.4,
    )
    sim.add_particles(1.0e-9, distr, npart)
    return sim


def linear_lattice():
    """Linear elements, incl. a bend, a solenoid, and a thin skew multipole."""
    return elements.KnownElementsList(
        [
            elements.Drift(name="d1", ds=0.5, nslice=2),
            elements.Quad(name="qf", ds=0.3, k=3.0, nslice=3),
            elements.Sbend(name="b1", ds=1.0, rc=4.0, nslice=4),
            elements.Multipole(name="m1", multipole=2, K_normal=0.5, K_skew=0.2),
            elements.Sol(name="s1", ds=0.4, ks=1.5, nslice=2),
            elements.Quad(name="qd", ds=0.3, k=-3.0, nslice=3),
            elements.Drift(name="d2", ds=0.5, nslice=2),
        ]
    )


def track(sim, lattice, stop=None):
    """Push the beam slice by slice through the lattice, up to element index stop"""
    for el in lattice[:stop]:
        sim.beam.ref.sedge = sim.beam.ref.s
        for _ in range(el.nslice):
            el.push(sim.beam)


def test_survey_geometry():
    """The reference trajectory through a sector bend"""
    sim = make_sim(0)
    ref = sim.beam.ref
    rc = 2.0
    angle = math.pi / 6
    lattice = elements.KnownElementsList(
        [
            elements.Drift(name="d1", ds=1.0),
            elements.Sbend(name="b1", ds=rc * angle, rc=rc, nslice=6),
            elements.Drift(name="d2", ds=1.0),
        ]
    )
    s_before = ref.s
    df = lattice.survey(ref)
    assert ref.s == s_before  # the reference particle is not changed
    assert len(df) == 1 + 1 + 6 + 1
    assert len(lattice.survey(ref, slices=False)) == 1 + 3

    end = df.iloc[-1]
    theta = end["theta"]
    assert abs(theta) == pytest.approx(angle, rel=1.0e-12)
    sign = math.copysign(1.0, theta)
    x = sign * rc * (1.0 - math.cos(angle)) + math.sin(theta)
    z = 1.0 + rc * math.sin(angle) + math.cos(theta)
    assert end["x"] == pytest.approx(x, rel=1.0e-12)
    assert end["z"] == pytest.approx(z, rel=1.0e-12)
    assert end["s"] == pytest.approx(2.0 + rc * angle, rel=1.0e-12)
    sim.finalize()


def test_survey_envelope_matches_tracking():
    """For linear elements, envelope tracking propagates the beam covariance exactly"""
    sim = make_sim()
    lattice = linear_lattice()
    df = lattice.survey(sim.beam.ref, covariance=sim.beam.covariance_matrix())

    track(sim, lattice)
    cm = sim.beam.covariance_matrix()
    rtol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-12
    for j, name in enumerate(SIGMAS, start=1):
        assert df[name].iloc[-1] == pytest.approx(math.sqrt(cm[j, j]), rel=rtol)
    # the beam sizes change along the lattice
    assert df["sigma_x"].max() > 2.0 * df["sigma_x"].min()
    sim.finalize()


@pytest.mark.parametrize("where", ["middle", "end"])
def test_survey_anchored(where):
    """A survey from a reference particle inside or at the end of the lattice"""
    sim = make_sim()
    lattice = linear_lattice()
    start = lattice.survey(sim.beam.ref, covariance=sim.beam.covariance_matrix())

    # middle: after the bend, before the thin multipole at the same position
    stop, element = (3, lattice[2]) if where == "middle" else (None, None)
    track(sim, lattice, stop)
    other = lattice.survey(
        sim.beam.ref, covariance=sim.beam.covariance_matrix(), element=element
    )

    rtol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-10
    for col in ("s", "x", "y", "z", "theta", "pt", *SIGMAS):
        expected = start[col].to_numpy()
        scale = np.max(np.abs(expected))
        if scale == 0.0:  # e.g., y in a horizontal lattice
            np.testing.assert_array_equal(other[col].to_numpy(), expected)
        else:
            np.testing.assert_allclose(
                other[col].to_numpy() / scale, expected / scale, rtol=0.0, atol=rtol
            )
    sim.finalize()


def test_survey_ambiguous_position():
    """Inside the lattice at zero-length elements that change the beam, the element is needed"""
    sim = make_sim()
    lattice = linear_lattice()
    track(sim, lattice, 3)  # after the bend, before the thin multipole
    cm = sim.beam.covariance_matrix()
    with pytest.raises(ValueError, match="m1"):
        lattice.survey(sim.beam.ref, covariance=cm)
    lattice.survey(sim.beam.ref, covariance=cm, element=2)

    # a marker at the same position does not change anything
    markers = elements.KnownElementsList(
        [
            elements.Drift(name="d1", ds=0.5, nslice=2),
            elements.Marker(name="mark"),
            elements.Drift(name="d2", ds=0.5, nslice=2),
        ]
    )
    sim2 = make_sim()
    track(sim2, markers, 1)
    markers.survey(sim2.beam.ref, covariance=sim2.beam.covariance_matrix())
    sim2.finalize()
    sim.finalize()


def test_survey_errors():
    """A reference particle between slice boundaries, and elements without a linear map"""
    sim = make_sim()
    lattice = linear_lattice()
    ref = sim.beam.ref.copy()
    ref.s = 0.1  # inside the first slice
    with pytest.raises(ValueError, match="slice boundary"):
        lattice.survey(ref)
    with pytest.raises(ValueError, match="not in this lattice"):
        lattice.survey(sim.beam.ref, element=elements.Drift(ds=1.0))

    programmable = elements.KnownElementsList(
        [elements.Programmable(name="prog", ds=0.5)]
    )
    with pytest.raises(RuntimeError, match="prog"):
        programmable.survey(sim.beam.ref, covariance=sim.beam.covariance_matrix())
    sim.finalize()


def bends_lattice(cells):
    """Drifts and 15 degree bends, in slices of 0.1 m"""
    angle = math.radians(15.0)
    lattice = elements.KnownElementsList()
    for _ in range(cells):
        lattice.extend(
            [
                elements.Drift(name="d", ds=1.0, nslice=10),
                elements.Sbend(name="b", ds=1.0, rc=1.0 / angle, nslice=10),
            ]
        )
    return lattice


def heading(t):
    """Angle of a direction in the horizontal plane [rad]"""
    return math.atan2(t[0], t[2])


def test_survey_smoothed_axes():
    """A camera that follows the smoothed orbit turns at a steady rate through bends"""
    from impactx.plot.Survey3D import frames, smoothed_axes

    sim = make_sim(0)
    df = bends_lattice(12).survey(sim.beam.ref)
    s = df["s"].to_numpy()
    _, t, n, v = frames(df)

    # without smoothing: the axes of the reference particle
    for row in (0, 15, len(df) - 1):
        axes = smoothed_axes(df, s[row], 0.0)
        for axis, local in zip(axes, (t[row], n[row], v[row])):
            assert axis == pytest.approx(local, abs=1.0e-12)

    sigma = 1.5
    s0 = np.linspace(s[0], s[-1], 241)
    smooth, local = [], []
    for position in s0:
        tc, nc, vc = smoothed_axes(df, position, sigma)
        assert np.linalg.norm(tc) == pytest.approx(1.0, rel=1.0e-12)
        assert nc[1] == pytest.approx(0.0, abs=1.0e-12)  # no roll
        assert np.cross(tc, nc) == pytest.approx(vc, abs=1.0e-12)
        smooth.append(heading(tc))
        local.append(heading(smoothed_axes(df, position, 0.0)[0]))
    turn, local_turn = np.diff(smooth), np.diff(local)

    # 12 cells of 15 degrees in 24 m; the bends are 2 m, about a sigma, apart
    rate = math.radians(15.0) / 2.0 * (s0[1] - s0[0])
    # far enough from the ends that the straight orbit outside has no weight
    inside = (s0[1:] > s[0] + 4.0 * sigma) & (s0[:-1] < s[-1] - 4.0 * sigma)
    assert np.abs(turn[inside]) == pytest.approx(rate, rel=1.0e-3)
    # the orbit continues straight before and after the lattice: the camera turns
    # slower there, without a jump
    assert np.max(np.abs(turn)) <= rate * (1.0 + 1.0e-3)
    assert np.max(np.abs(np.diff(turn))) < 0.05 * np.max(np.abs(np.diff(local_turn)))

    # a straight lattice
    straight = elements.KnownElementsList([elements.Drift(ds=2.0, nslice=20)])
    df_straight = straight.survey(sim.beam.ref)
    _, t_straight, _, _ = frames(df_straight)
    assert smoothed_axes(df_straight, 1.0, sigma)[0] == pytest.approx(
        t_straight[0], abs=1.0e-12
    )
    sim.finalize()


def test_survey_smoothed_axes_ring():
    """In a ring, the smoothed orbit wraps around: no jump from turn to turn"""
    from impactx.plot.Survey3D import smoothed_axes

    sim = make_sim(0)
    df = bends_lattice(24).survey(sim.beam.ref)  # 360 degrees
    s = df["s"].to_numpy()
    sigma = 1.5
    start, end = (smoothed_axes(df, position, sigma) for position in (s[0], s[-1]))
    for a, b in zip(start, end):
        assert a == pytest.approx(b, abs=1.0e-9)

    # the camera turns at the rate inside the ring, also at the start of a turn
    ds = 0.1
    rate = math.radians(15.0) / 2.0 * ds
    turn = heading(smoothed_axes(df, s[0] + ds, sigma)[0]) - heading(start[0])
    assert abs(turn) == pytest.approx(rate, rel=1.0e-3)
    sim.finalize()


def test_survey_orbit_position():
    """Positions on the orbit: interpolated, straight past the ends, around a ring"""
    from impactx.plot.Survey3D import frames, orbit_position

    sim = make_sim(0)
    df = bends_lattice(1).survey(sim.beam.ref)  # 1 m drift, then a 15 degree bend
    p, t, _, _ = frames(df)
    s = df["s"].to_numpy()

    # straight before and after the lattice
    assert orbit_position(df, s[0] - 0.5) == pytest.approx(p[0] - 0.5 * t[0], abs=1e-12)
    assert orbit_position(df, s[-1] + 0.5) == pytest.approx(
        p[-1] + 0.5 * t[-1], abs=1e-12
    )
    # in the bend, between two survey rows 0.1 m apart: on the arc, within the
    # sagitta of a chord of 0.1 m
    rc = 1.0 / math.radians(15.0)
    phi = 0.55 / rc
    sign = math.copysign(1.0, p[-1][0])
    expected = [sign * rc * (1.0 - math.cos(phi)), 0.0, 1.0 + rc * math.sin(phi)]
    assert orbit_position(df, 1.55) == pytest.approx(expected, abs=0.1**2 / (8 * rc))

    # a ring wraps around
    ring = bends_lattice(24).survey(sim.beam.ref)  # 360 degrees
    length = ring["s"].iloc[-1] - ring["s"].iloc[0]
    for position in (0.3, 7.25):
        assert orbit_position(ring, position + length) == pytest.approx(
            orbit_position(ring, position), abs=1e-9
        )
        assert orbit_position(ring, position - length) == pytest.approx(
            orbit_position(ring, position), abs=1e-9
        )
    sim.finalize()


def pipes_lattice():
    """A drift, a quad, a drift, a bend, a drift"""
    return elements.KnownElementsList(
        [
            elements.Drift(name="d1", ds=1.0, nslice=4),
            elements.Quad(name="q1", ds=0.3, k=2.0, nslice=3),
            elements.Drift(name="d2", ds=0.5, nslice=2),
            elements.Sbend(name="b1", ds=1.0, rc=4.0, nslice=4),
            elements.Drift(name="d3", ds=0.5, nslice=2),
        ]
    )


def test_survey_geometry_pipes():
    """Cut open lattices with a bore draw their drifts as beam pipes"""
    from impactx.plot.Survey3DGeometry import Geometry

    sim = make_sim(0)
    lattice = pipes_lattice()
    df = lattice.survey(sim.beam.ref)

    def geometry(**kwargs):
        return Geometry(lattice, df, None, 0.0, 1.0, 0.1, "cern-lhc", **kwargs)

    label = Geometry.PIPE_LABEL
    cut = geometry(bore="round", cutaway=True)
    assert label in cut.legend
    assert cut.legend[label] == Geometry.PIPE_COLOR
    # the drifts add surfaces to the quad and the bend
    no_drifts = elements.KnownElementsList([el for el in lattice if el.name[0] != "d"])
    without = Geometry(
        no_drifts,
        no_drifts.survey(sim.beam.ref),
        None,
        0.0,
        1.0,
        0.1,
        "cern-lhc",
        bore="round",
        cutaway=True,
    )
    assert len(cut._coords) > len(without._coords)
    # the inside of the pipes is the bore of the elements
    pipe = cut._kinds == cut.colors.index(Geometry.PIPE_COLOR) + 1
    p = df[["x", "y", "z"]].to_numpy()
    distance = np.min(
        np.linalg.norm(cut._coords[pipe][:, None, :] - p[None, :, :], axis=2), axis=1
    )
    radius = Geometry.BORE * 0.1
    assert distance.min() == pytest.approx(radius, rel=0.05)
    assert distance.max() <= radius + Geometry.PIPE_WALL * 0.1 + 1.0e-9

    # no pipes: translucent elements, or solid cut open ones without a bore
    assert label not in geometry(bore="round", cutaway=False).legend
    assert label not in geometry(bore=None, cutaway=True).legend
    sim.finalize()


def test_survey_geometry_opaque():
    """Opaque elements are closed around the bore, with drifts as beam pipes: a tunnel"""
    from impactx.plot.Survey3D import SURVEY_TOPOLOGY, frames
    from impactx.plot.Survey3DGeometry import Geometry

    sim = make_sim(0)
    lattice = pipes_lattice()
    # a zero-length element between the second drift and the bend
    lattice.insert(
        3, elements.Multipole(name="m1", multipole=2, K_normal=0.0, K_skew=0.0)
    )
    df = lattice.survey(sim.beam.ref)
    a = 0.1

    def geometry(**kwargs):
        return Geometry(lattice, df, None, 0.0, 1.0, a, "cern-lhc", **kwargs)

    p, t, _, v = frames(df)

    def local(points):
        """Distance from the orbit and height above it, transverse to the nearest
        orbit point"""
        nearest = np.argmin(
            np.linalg.norm(points[:, None, :] - p[None, :, :], axis=2), axis=1
        )
        offset = points - p[nearest]
        along = np.sum(offset * t[nearest], axis=1)
        offset = offset - along[:, None] * t[nearest]
        return np.linalg.norm(offset, axis=1), np.sum(offset * v[nearest], axis=1)

    label = Geometry.PIPE_LABEL
    tunnel = geometry(bore="round", opaque=True)
    cut = geometry(bore="round", cutaway=True)
    for g in (tunnel, cut):
        assert label in g.legend
    shape = tunnel.blueprint()["topologies"][SURVEY_TOPOLOGY]["elements"]["shape"]
    assert shape == "quad"

    radius = Geometry.BORE * a
    pipe = tunnel._kinds == tunnel.colors.index(Geometry.PIPE_COLOR) + 1
    distance, height = local(tunnel._coords[pipe])
    assert distance.min() == pytest.approx(radius, rel=0.05)
    assert distance.max() <= radius + Geometry.PIPE_WALL * a + 1.0e-9
    # closed: the pipes go all around the bore, cut open only below the axis
    assert height.max() >= 0.99 * radius
    assert height.min() <= -0.99 * radius
    _, cut_height = local(cut._coords)
    # the cut of the zero-length element is a little above the others
    assert cut_height.max() <= Geometry.SLAB_LIFT * a + 1.0e-9
    _, tunnel_height = local(tunnel._coords[~pipe])
    assert tunnel_height.max() >= 0.99 * radius

    # the zero-length element is a collar, inside the bore of its neighbors
    multipole = tunnel._kinds == tunnel.colors.index(tunnel.legend["Multipole"]) + 1
    distance, _ = local(tunnel._coords[multipole])
    assert distance.min() == pytest.approx(radius - Geometry.COLLAR * a, rel=0.05)
    sim.finalize()


def test_survey_repeated_element():
    """An element that is in the lattice more than once is the one at the position of
    the reference particle, e.g., sim.tracking_element in a tracking hook"""
    sim = make_sim(0)
    drift = elements.Drift(name="d", ds=0.5, nslice=2)
    bend = elements.Sbend(name="b", ds=1.0, rc=4.0, nslice=4)

    def push_to(lattice, stop):
        """The reference particle after the elements before index stop"""
        ref = sim.beam.ref.copy()
        for el in list(lattice)[:stop]:
            ref.sedge = ref.s
            for _ in range(el.nslice):
                push(ref, el)
        return ref

    # the same drift at the indices 0, 2 and 4, around a marker
    lattice = elements.KnownElementsList(
        [drift, bend, drift, elements.Marker(name="m"), drift, bend]
    )
    for stop, index in ((1, 0), (2, 2), (4, 4), (5, 4)):
        ref = push_to(lattice, stop)
        pd.testing.assert_frame_equal(
            lattice.survey(ref, element=drift), lattice.survey(ref, element=index)
        )
    # around a thin quadrupole, before or after it is ambiguous: the index is needed
    multipole = elements.Multipole(name="m", multipole=2, K_normal=0.5, K_skew=0.0)
    focused = elements.KnownElementsList([drift, bend, drift, multipole, drift])
    ref = push_to(focused, 3)
    with pytest.raises(ValueError, match="more than once"):
        focused.survey(ref, element=drift)
    for index in (2, 4):
        focused.survey(ref, element=index)
    sim.finalize()


def test_survey_overlapping_slabs():
    """Cut open zero-length elements that overlap have no coinciding surfaces"""
    from impactx.plot.Survey3DGeometry import Geometry

    sim = make_sim(0)
    lattice = elements.KnownElementsList(
        [
            elements.Drift(name="d1", ds=0.5, nslice=2),
            elements.Multipole(name="m", multipole=2, K_normal=0.0, K_skew=0.0),
            elements.Kicker(name="k", xkick=0.0, ykick=0.0),
            elements.Drift(name="d2", ds=0.5, nslice=2),
        ]
    )
    a = 0.1
    df = lattice.survey(sim.beam.ref)
    g = Geometry(lattice, df, None, 0.0, 1.0, a, "cern-lhc", cutaway=True)

    def extent(kind):
        """Thickness along the straight orbit, full width, and height of the cut"""
        points = g._coords[g._kinds == g.colors.index(g.legend[kind]) + 1]
        size = points.max(axis=0) - points.min(axis=0)
        return size[2], size[0], points[:, 1].max()

    thin, wide, cut = (Geometry.THIN_OPAQUE, Geometry.SLAB_WIDER, Geometry.SLAB_LIFT)
    half_width = 0.9 * a
    expected = {
        "Multipole": (thin * a, 2 * half_width * (1 + wide), cut * a),
        "Kicker": (
            thin * (1 - Geometry.SLAB_THINNER) * a,
            2 * half_width * (1 + 2 * wide),
            2 * cut * a,
        ),
    }
    for kind, values in expected.items():
        assert extent(kind) == pytest.approx(values, rel=1e-9), kind
    sim.finalize()


def test_survey_geometry_cache():
    """Elements reuse their meshes from a cache, unless their position, size or bore
    change, or the settings of the geometry"""
    from impactx.plot.Survey3DGeometry import Geometry

    sim = make_sim()
    lattice = pipes_lattice()
    start = lattice.survey(sim.beam.ref, covariance=sim.beam.covariance_matrix())
    # the same lattice, surveyed from the reference particle after its second element
    for el in lattice[:2]:
        sim.beam.ref.sedge = sim.beam.ref.s
        for _ in range(el.nslice):
            el.push(sim.beam)
    later = lattice.survey(
        sim.beam.ref, covariance=sim.beam.covariance_matrix(), element=1
    )
    a = 0.1

    def geometry(df, magnify=1.0, size=a, cache=None):
        return Geometry(
            lattice, df, df, 1.0, magnify, size, "cern-lhc", cutaway=True, cache=cache
        )

    def meshes(cache):
        return {k: v[3] for k, v in cache.items() if isinstance(k, int)}

    # the envelope fits into the smallest bore: anchored elsewhere, the geometry is
    # the same within roundoff, and all elements and pipes reuse their meshes
    cache = {}
    geometry(start, cache=cache)
    before = meshes(cache)
    assert len(before) == len(lattice)
    reused = geometry(later, cache=cache)
    assert all(meshes(cache)[k] is before[k] for k in before)
    fresh = geometry(later)
    np.testing.assert_allclose(reused._coords, fresh._coords, rtol=0.0, atol=1e-9 * a)
    np.testing.assert_array_equal(reused._conn, fresh._conn)
    np.testing.assert_array_equal(reused._kinds, fresh._kinds)

    # a larger envelope widens the bores: all elements are rebuilt
    wide = geometry(later, magnify=2.0e3, cache=cache)
    assert all(meshes(cache)[k] is not before[k] for k in before)
    np.testing.assert_array_equal(wide._coords, geometry(later, magnify=2.0e3)._coords)

    # an element size that differs by roundoff, e.g., relative to the extent of the
    # lattice, keeps the cache; another one, or other settings, start over
    before = meshes(cache)
    geometry(later, magnify=2.0e3, size=a * (1.0 + 1.0e-14), cache=cache)
    assert all(meshes(cache)[k] is before[k] for k in before)
    geometry(later, size=2.0 * a, cache=cache)
    assert cache["element_size"] == 2.0 * a
    assert all(meshes(cache)[k] is not before[k] for k in before)
    before = meshes(cache)
    Geometry(
        lattice, later, later, 1.0, 1.0, 2.0 * a, "cern-lhc", opaque=True, cache=cache
    )
    assert all(meshes(cache)[k] is not before[k] for k in before)
    sim.finalize()


def test_survey_wires_camera():
    """Wires close to a camera in the lattice are thinner, in proportion to the distance"""
    from impactx.plot.Survey3DGeometry import Geometry

    sim = make_sim()
    lattice = pipes_lattice()
    df = lattice.survey(sim.beam.ref, covariance=sim.beam.covariance_matrix())
    geometry = Geometry(lattice, df, df, 1.0, 100.0, 0.1, "cern-lhc")
    p = geometry._p

    def orbit_half_thickness():
        """Half thickness of the orbit wire at each orbit point: the first wire"""
        corners = geometry._wires[0][: 4 * len(p)].reshape(len(p), 4, 3)
        return np.linalg.norm(corners - p[:, None, :], axis=2).max(axis=1) / math.sqrt(
            2
        )

    radius = 0.002
    geometry.add_wires(radius, 0.3)
    assert orbit_half_thickness() == pytest.approx(radius, rel=1e-9)
    num_cells = len(geometry._wires[1])

    # a camera inside the lattice, 0.03 m above the middle of the orbit
    eye = p[len(p) // 2] + 0.03 * geometry._v[len(p) // 2]
    geometry.add_wires(radius, 0.3, camera=(eye, 0.5))
    distance = np.linalg.norm(p - eye[None, :], axis=1)
    assert orbit_half_thickness() == pytest.approx(radius * distance / 0.5, rel=1e-9)
    # the whole cage, in the lab frame: wherever the camera is
    assert len(geometry._wires[1]) == num_cells
    sim.finalize()


def racetrack_lattice(gap=0.0, curvature=1.0):
    """A ring of two 10 m straights and two arcs of 12 bends of 15 degrees

    The second straight is longer by gap. The sign of curvature sets the sense of rotation.
    """
    rc = curvature / math.radians(15.0)
    lattice = elements.KnownElementsList()
    for length in (10.0, 10.0 + gap):
        lattice.append(elements.Drift(name="straight", ds=length, nslice=4))
        lattice.extend(
            elements.Sbend(name="b", ds=1.0, rc=rc, nslice=2) for _ in range(12)
        )
    return lattice


# the sign of the closing error of a ring must not matter
@pytest.mark.parametrize("gap", [-1.0e-9, 1.0e-9])
def test_survey_overview_camera_ring(gap):
    """The overview of a ring stays in place while the beam goes around"""
    from impactx.extensions.InSituSurvey3D import _Survey3D

    sim = make_sim(0)
    lattice = racetrack_lattice(gap)
    ref = sim.beam.ref.copy()

    def camera():
        actions = _Survey3D(lattice, ref=ref, envelope=0.0).actions()
        scenes = next(a["scenes"] for a in actions if a["action"] == "add_scenes")
        return scenes["survey3d"]["renders"]["r1"]["camera"]

    tolerance = 1.0e-3 if Config.precision == "SINGLE" else 1.0e-6
    first = camera()
    for element in lattice[:-1]:
        for _ in range(element.nslice):
            push(ref, element)
        view = camera()
        for key in ("position", "look_at", "up"):
            assert view[key] == pytest.approx(first[key], abs=tolerance)
    sim.finalize()


@pytest.mark.parametrize("camera", ["beam", "racecar"])
def test_survey_camera_inside_bounds(camera):
    """Close to the beam, the camera is inside the bounds of the ray traced surfaces and
    wires, so that Viskores traces all pixels, also with corners of the bounds behind
    the camera"""
    from impactx.extensions.InSituSurvey3D import _Survey3D
    from impactx.plot.Survey3D import SURVEY_TOPOLOGY, WIRES_TOPOLOGY

    sim = make_sim(0)
    vis = _Survey3D(
        racetrack_lattice(),
        ref=sim.beam.ref,
        envelope=0.0,
        camera=camera,
        element_size=0.1,
        cutaway=True,
    )
    actions = vis.actions()
    scenes = next(a["scenes"] for a in actions if a["action"] == "add_scenes")
    position = scenes["survey3d"]["renders"]["r1"]["camera"]["position"]
    mesh = vis.geometry()
    for topology in (SURVEY_TOPOLOGY, WIRES_TOPOLOGY):
        coords = mesh["coordsets"][f"{topology}_coords"]["values"]
        for axis, p in zip("xyz", position):
            assert coords[axis].min() <= p <= coords[axis].max()
    sim.finalize()


@pytest.mark.parametrize("curvature", [1.0, -1.0])
def test_survey_overview_axis(curvature):
    """The direction to the right in the overview of a beamline and of rings"""
    from impactx.plot.Survey3D import frames, overview_axis

    sim = make_sim(0)

    # a beamline: along its net travel
    df = bends_lattice(3).survey(sim.beam.ref)
    p, _, _, _ = frames(df)
    axis = overview_axis(df)
    assert np.linalg.norm(axis) == pytest.approx(1.0)
    assert axis @ (p[-1] - p[0])[[0, 2]] > 0.0

    # a ring: along its long axis, with the motion at its start
    axis = overview_axis(racetrack_lattice(curvature=curvature).survey(sim.beam.ref))
    assert axis == pytest.approx([0.0, 1.0], abs=0.01)  # the first straight is along z
    sim.finalize()

    # a circle has no long axis: the direction of motion at its start
    angle = curvature * np.linspace(0.0, 2.0 * np.pi, 10001)
    circle = pd.DataFrame(
        {
            "x": np.cos(angle),
            "y": 0.0,
            "z": np.sin(angle),
            "px": -np.sin(angle),
            "py": 0.0,
            "pz": np.cos(angle),
        }
    )
    assert overview_axis(circle) == pytest.approx([0.0, 1.0])
