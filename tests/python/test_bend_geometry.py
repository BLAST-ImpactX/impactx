#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""The geometry of a bend: ds plus one of rc, phi or B, or phi together with B."""

import math
import warnings

import numpy as np
import pytest

from impactx import ImpactX, RefPart, elements

BENDS = {
    "Sbend": lambda **geometry: elements.Sbend(**geometry),
    "ExactSbend": lambda **geometry: elements.ExactSbend(**geometry),
    "CFbend": lambda **geometry: elements.CFbend(k=0.3, **geometry),
}


def electron(kin_energy_MeV=1.0e3):
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)
    return ref


def geometries(ds, rc, ref):
    """The same bend of length ds and radius rc, specified each way"""
    phi = math.degrees(ds / rc)
    B = ref.rigidity_Tm / rc
    return {
        "rc": dict(ds=ds, rc=rc),
        "phi": dict(ds=ds, phi=phi),
        "B": dict(ds=ds, B=B),
        # a fixed angle: ds does not set the arc length, so pick a different one
        "phi_B": dict(ds=2.0 * ds, phi=phi, B=B),
    }


def track_reference(element, kin_energy_MeV=1.0e3):
    """The reference particle after the element"""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.init_grids()
    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)
    sim.lattice.append(element)
    sim.track_reference(ref)
    result = dict(x=ref.x, z=ref.z, px=ref.px, pz=ref.pz, t=ref.t)
    sim.finalize()
    return result


@pytest.mark.parametrize("kind", BENDS.keys())
@pytest.mark.parametrize("rc", [10.0, -10.0])
def test_every_specification_gives_the_same_bend(kind, rc):
    ref = electron()
    maps = {}
    orbits = {}
    for name, geometry in geometries(0.5, rc, ref).items():
        bend = BENDS[kind](**geometry)
        assert bend.signed_rc(ref) == pytest.approx(rc, rel=1e-12), name
        maps[name] = np.array(bend.transfer_map(ref))
        orbits[name] = track_reference(bend)

    for name in ["phi", "B", "phi_B"]:
        np.testing.assert_allclose(maps[name], maps["rc"], rtol=1e-10, atol=1e-14)
        # phi_B has its own ds: s differs, but not the orbit
        for key, value in orbits["rc"].items():
            assert orbits[name][key] == pytest.approx(value, rel=1e-10, abs=1e-14), (
                name,
                key,
            )


def test_a_fixed_angle_is_kept_at_any_energy():
    """A cyclotron half turn: the radius grows with the energy, the angle stays 180 degrees"""
    radii = []
    for kin_energy_MeV in [1.0, 10.0]:
        half_turn = elements.ExactSbend(ds=0.25, phi=180.0, B=1.0)
        orbit = track_reference(half_turn, kin_energy_MeV)
        pz0 = electron(kin_energy_MeV).beta_gamma
        assert orbit["pz"] == pytest.approx(-pz0, rel=1e-12)
        assert orbit["px"] == pytest.approx(0.0, abs=1e-12 * pz0)
        radii.append(abs(half_turn.signed_rc(electron(kin_energy_MeV))))
    assert radii[1] > radii[0]


def test_a_field_bends_by_the_charge():
    electron_ref = electron()
    positron_ref = RefPart()
    positron_ref.set_species("positron").set_kin_energy_MeV(1.0e3)

    by_field = elements.Sbend(ds=0.5, B=0.5)
    assert by_field.signed_rc(positron_ref) == pytest.approx(
        -by_field.signed_rc(electron_ref), rel=1e-12
    )

    by_radius = elements.Sbend(ds=0.5, rc=10.0)
    assert by_radius.signed_rc(positron_ref) == by_radius.signed_rc(electron_ref)


@pytest.mark.parametrize("kind", BENDS.keys())
@pytest.mark.parametrize(
    "geometry",
    [
        dict(),
        dict(rc=10.0, phi=2.0),
        dict(rc=10.0, B=0.5),
        dict(rc=10.0, phi=2.0, B=0.5),
    ],
    ids=["none", "rc_phi", "rc_B", "rc_phi_B"],
)
def test_one_specification_is_required(kind, geometry):
    with pytest.raises(ValueError, match="rc"):
        BENDS[kind](ds=0.5, **geometry)


@pytest.mark.parametrize("kind", BENDS.keys())
def test_legacy_zero_field_with_an_angle_is_the_angle_alone(kind):
    with pytest.warns(DeprecationWarning, match="B=0.0 together with phi"):
        legacy = BENDS[kind](ds=0.5, phi=3.0, B=0.0)
    assert legacy.B is None
    assert legacy.to_dict(**_in_degrees(kind)) == BENDS[kind](ds=0.5, phi=3.0).to_dict(
        **_in_degrees(kind)
    )


def _in_degrees(kind):
    """ExactSbend.to_dict() gives phi in radians unless asked for degrees"""
    return dict(in_degrees=True) if kind == "ExactSbend" else {}


@pytest.mark.parametrize("kind", ["Sbend", "CFbend"])
def test_parameters_that_do_not_specify_the_bend_are_none(kind):
    bend = BENDS[kind](ds=0.5, phi=3.0, B=0.5)
    assert (bend.rc, bend.phi, bend.B) == (None, 3.0, 0.5)
    d = bend.to_dict()
    assert (d["rc"], d["phi"], d["B"]) == (None, 3.0, 0.5)


@pytest.mark.parametrize("kind", BENDS.keys())
def test_a_property_changes_a_parameter_of_the_specification(kind):
    bend = BENDS[kind](ds=0.5, rc=10.0)
    bend.rc = 12.0
    assert bend.signed_rc(electron()) == 12.0

    with pytest.raises(ValueError, match="set_geometry"):
        bend.B = 0.5
    with pytest.raises(ValueError, match="set_geometry"):
        bend.rc = None
    assert (bend.rc, bend.phi, bend.B) == (12.0, None, None)


@pytest.mark.parametrize("kind", BENDS.keys())
def test_set_geometry_changes_the_specification(kind):
    ref = electron()
    bend = BENDS[kind](ds=0.5, rc=10.0)

    bend.set_geometry(rc=None, B=0.5)
    assert (bend.rc, bend.B) == (None, 0.5)
    assert bend.signed_rc(ref) == pytest.approx(ref.rigidity_Tm / 0.5, rel=1e-12)

    # the remaining parameters are kept: adding phi gives a fixed angle
    bend.set_geometry(phi=90.0)
    assert bend.B == 0.5
    assert bend.to_dict(**_in_degrees(kind))["phi"] == pytest.approx(90.0)

    # a rejected change leaves the bend as it was
    with pytest.raises(ValueError):
        bend.set_geometry(rc=10.0)
    assert (bend.rc, bend.B) == (None, 0.5)

    with pytest.raises(TypeError, match="unexpected keyword"):
        bend.set_geometry(angle=10.0)


@pytest.mark.parametrize("kind", BENDS.keys())
def test_copy_changes_the_specification_together(kind):
    bend = BENDS[kind](ds=0.5, rc=10.0, name="b")
    by_angle = bend.copy(rc=None, phi=3.0)
    assert by_angle.to_dict(**_in_degrees(kind))["phi"] == pytest.approx(3.0)
    assert by_angle.rc is None
    assert bend.rc == 10.0  # the original is unchanged

    with pytest.raises(ValueError):
        bend.copy(phi=3.0)  # rc and phi


@pytest.mark.parametrize("kind", BENDS.keys())
@pytest.mark.parametrize("spec", ["rc", "phi", "B", "phi_B"])
def test_dicts_round_trip(kind, spec):
    geometry = geometries(0.5, 10.0, electron())[spec]
    lattice = elements.KnownElementsList()
    lattice.append(BENDS[kind](name="b", **geometry))

    restored = elements.KnownElementsList()
    restored.from_dicts(lattice.to_dicts())
    assert restored[0].to_dict(**_in_degrees(kind)) == pytest.approx(
        lattice[0].to_dict(**_in_degrees(kind))
    )


@pytest.mark.parametrize("kind", BENDS.keys())
@pytest.mark.parametrize(
    "geometry", [dict(rc=0.0), dict(phi=0.0), dict(B=0.0)], ids=["rc", "phi", "B"]
)
def test_a_straight_bend_is_a_drift_of_its_orbit(kind, geometry):
    ref = electron()
    straight = BENDS[kind](ds=0.5, **geometry)
    assert math.isinf(straight.signed_rc(ref))

    orbit = track_reference(straight)
    drift = track_reference(elements.Drift(ds=0.5))
    for key, value in drift.items():
        assert orbit[key] == pytest.approx(value, rel=1e-12, abs=1e-14), key


def test_a_straight_cfbend_is_a_quadrupole():
    ref = electron()
    np.testing.assert_allclose(
        np.array(elements.CFbend(ds=0.5, rc=0.0, k=0.3).transfer_map(ref)),
        np.array(elements.Quad(ds=0.5, k=0.3).transfer_map(ref)),
        rtol=1e-12,
        atol=1e-15,
    )


def test_input_file_geometry(tmp_path):
    """Input files take the same rc, phi and B keys"""
    inputs = tmp_path / "inputs"
    inputs.write_text(
        "lattice.elements = geo_sb_rc geo_sb_phi geo_esb_b geo_esb_phib geo_cf_phi"
        " geo_td_b geo_de_b\n"
        "geo_sb_rc.type = sbend\ngeo_sb_rc.ds = 0.5\ngeo_sb_rc.rc = 10.0\n"
        "geo_sb_phi.type = sbend\ngeo_sb_phi.ds = 0.5\ngeo_sb_phi.phi = 3.0\n"
        "geo_esb_b.type = sbend_exact\ngeo_esb_b.ds = 0.5\ngeo_esb_b.B = 0.5\n"
        "geo_esb_phib.type = sbend_exact\ngeo_esb_phib.ds = 0.25\n"
        "geo_esb_phib.phi = 180.0\ngeo_esb_phib.B = 1.0\n"
        "geo_cf_phi.type = cfbend\ngeo_cf_phi.ds = 0.5\ngeo_cf_phi.phi = 3.0\n"
        "geo_cf_phi.k = 0.3\n"
        "geo_td_b.type = thin_dipole\ngeo_td_b.theta = 1.0\ngeo_td_b.B = 0.5\n"
        "geo_de_b.type = dipedge\ngeo_de_b.psi = 0.1\ngeo_de_b.g = 0.05\n"
        "geo_de_b.B = 0.5\n"
    )

    sim = ImpactX()
    sim.load_inputs_file(str(inputs))
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.init_grids()
    sim.init_lattice_elements_from_inputs()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # ExactSbend.to_dict()
        specs = [(el.rc, el.to_dict().get("B"), el.name) for el in sim.lattice]
    assert specs == [
        (10.0, None, "geo_sb_rc"),
        (None, None, "geo_sb_phi"),
        (None, 0.5, "geo_esb_b"),
        (None, 1.0, "geo_esb_phib"),
        (None, None, "geo_cf_phi"),
        (None, 0.5, "geo_td_b"),
        (None, 0.5, "geo_de_b"),
    ]
    assert sim.lattice[1].phi == pytest.approx(3.0)
    assert sim.lattice[4].phi == pytest.approx(3.0)

    sim.finalize()


# thin bends and edges ########################################################

THIN_AND_EDGE = {
    "ThinDipole": lambda **radius: elements.ThinDipole(theta=1.0, **radius),
    "DipEdge": lambda **radius: elements.DipEdge(psi=0.1, g=0.05, **radius),
}


@pytest.mark.parametrize("kind", THIN_AND_EDGE.keys())
@pytest.mark.parametrize("rc", [10.0, -10.0])
def test_thin_and_edge_radius_or_field(kind, rc):
    ref = electron()
    by_radius = THIN_AND_EDGE[kind](rc=rc)
    by_field = THIN_AND_EDGE[kind](B=ref.rigidity_Tm / rc)
    assert (by_field.rc, by_radius.B) == (None, None)
    np.testing.assert_allclose(
        np.array(by_field.transfer_map(ref)),
        np.array(by_radius.transfer_map(ref)),
        rtol=1e-12,
        atol=1e-15,
    )


@pytest.mark.parametrize(
    "construct",
    [
        lambda: elements.ThinDipole(theta=1.0),
        lambda: elements.ThinDipole(theta=1.0, rc=10.0, B=0.5),
        lambda: elements.DipEdge(psi=0.1, g=0.05),
        lambda: elements.DipEdge(psi=0.1, g=0.05, rc=10.0, B=0.5),
    ],
    ids=["thin_theta", "thin_rc_B", "edge_none", "edge_rc_B"],
)
def test_thin_and_edge_need_one_radius(construct):
    with pytest.raises(ValueError, match="rc"):
        construct()


def test_the_angle_of_a_thin_bend_is_theta():
    thin = elements.ThinDipole(theta=1.0, rc=10.0)
    thin.set_geometry(rc=None, B=0.5)
    assert (thin.rc, thin.B) == (None, 0.5)
    assert thin.to_dict(in_degrees=True)["theta"] == pytest.approx(1.0)
    with pytest.raises(TypeError, match="unexpected keyword"):
        thin.set_geometry(phi=2.0)
    with pytest.raises(ValueError):
        thin.set_geometry(theta=None)

    edge = elements.DipEdge(psi=0.1, g=0.05, rc=10.0)
    assert not hasattr(edge, "phi") and not hasattr(edge, "theta")
    with pytest.raises(TypeError, match="unexpected keyword"):
        edge.set_geometry(phi=2.0)


@pytest.mark.parametrize("kind", THIN_AND_EDGE.keys())
@pytest.mark.parametrize("radius", [dict(rc=0.0), dict(B=0.0)], ids=["rc", "B"])
def test_a_straight_thin_bend_or_edge_does_not_kick(kind, radius):
    ref = electron()
    np.testing.assert_array_equal(
        np.array(THIN_AND_EDGE[kind](**radius).transfer_map(ref)), np.eye(6)
    )


def test_a_straight_thin_bend_does_not_bend_the_orbit():
    orbit = track_reference(elements.ThinDipole(theta=1.0, rc=0.0))
    assert (orbit["px"], orbit["x"]) == (0.0, 0.0)


@pytest.mark.parametrize("kind", THIN_AND_EDGE.keys())
@pytest.mark.parametrize("radius", [dict(rc=10.0), dict(B=0.5)], ids=["rc", "B"])
def test_thin_and_edge_dicts_round_trip(kind, radius):
    lattice = elements.KnownElementsList()
    lattice.append(THIN_AND_EDGE[kind](name="b", **radius))

    restored = elements.KnownElementsList()
    restored.from_dicts(lattice.to_dicts())
    extra = dict(in_degrees=True) if kind == "ThinDipole" else {}
    assert restored[0].to_dict(**extra) == pytest.approx(lattice[0].to_dict(**extra))


# exact combined-function bends ###############################################


def exact_cfbend(**kwargs):
    """A combined-function bend with a quadrupole and a sextupole component"""
    kwargs.setdefault("k_normal", [0.0, 0.3, 0.5])
    kwargs.setdefault("k_skew", [0.0, 0.0, 0.1])
    return elements.ExactCFbend(mapsteps=4, **kwargs)


def track_beam(element, kin_energy_MeV=1.0e3):
    """A few electrons with spin after the element, as rows of (x, px, y, py, t, pt, sx, sy, sz)"""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.spin = True
    sim.init_grids()
    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)
    qm_eev = -1.0 / 0.51099895000 / 1e6  # electron charge/mass in e / eV
    offsets = [0.0, 1.0e-3, -2.0e-3]
    sim.beam.add_n_particles(
        offsets,
        [0.0, -1.0e-3, 5.0e-4],
        [0.0, 1.0e-4, 0.0],
        [0.0, 2.0e-4, -1.0e-4],
        [1.0e-4, 0.0, 0.0],
        [0.0, 1.0e-3, -1.0e-3],
        qm_eev,
        1.0e-12,
        sx=[1.0, 0.0, 0.6],
        sy=[0.0, 1.0, 0.0],
        sz=[0.0, 0.0, 0.8],
    )
    sim.lattice.append(element)
    sim.track_particles()
    df = sim.beam.to_df(local=True).sort_index()
    columns = [
        "position_x",
        "momentum_x",
        "position_y",
        "momentum_y",
        "position_t",
        "momentum_t",
        "spin_x",
        "spin_y",
        "spin_z",
    ]
    result = df[columns].to_numpy()
    sim.finalize()
    return result


@pytest.mark.parametrize("spec", ["rc", "phi", "B", "phi_B"])
@pytest.mark.parametrize("rc", [10.0, -10.0])
def test_exact_cfbend_geometry_is_the_dipole_coefficient(spec, rc):
    ref = electron()
    by_coefficient = exact_cfbend(ds=0.5, k_normal=[1.0 / rc, 0.3, 0.5])
    by_geometry = exact_cfbend(**geometries(0.5, rc, ref)[spec])

    assert by_geometry.signed_rc(ref) == pytest.approx(rc, rel=1e-12)
    np.testing.assert_allclose(
        np.array(by_geometry.transfer_map(ref)),
        np.array(by_coefficient.transfer_map(ref)),
        rtol=1e-10,
        atol=1e-14,
    )
    orbit, expected = track_reference(by_geometry), track_reference(by_coefficient)
    for key, value in expected.items():
        assert orbit[key] == pytest.approx(value, rel=1e-10, abs=1e-14), key
    np.testing.assert_allclose(
        track_beam(by_geometry), track_beam(by_coefficient), rtol=1e-9, atol=1e-14
    )


def test_exact_cfbend_dipole_is_given_once():
    with pytest.raises(ValueError, match=r"k_normal\[0\] must be 0"):
        exact_cfbend(ds=0.5, k_normal=[0.1, 0.3, 0.5], rc=10.0)

    by_coefficient = exact_cfbend(ds=0.5, k_normal=[0.1, 0.3, 0.5])
    assert (by_coefficient.rc, by_coefficient.phi, by_coefficient.B) == (
        None,
        None,
        None,
    )
    with pytest.raises(ValueError, match=r"k_normal\[0\] must be 0"):
        by_coefficient.set_geometry(rc=10.0)

    by_geometry = exact_cfbend(ds=0.5, rc=10.0)
    with pytest.raises(ValueError, match=r"k_normal\[0\] must be 0"):
        by_geometry.set_coefficients([0.1, 0.3, 0.5], [0.0, 0.0, 0.1])
    assert by_geometry.k_normal == [0.0, 0.3, 0.5]

    # copy() switches between both, in either direction
    switched = by_geometry.copy(
        rc=None, k_normal=[0.1, 0.3, 0.5], k_skew=[0.0, 0.0, 0.1]
    )
    assert (switched.rc, switched.k_normal[0]) == (None, 0.1)
    back = switched.copy(rc=10.0, k_normal=[0.0, 0.3, 0.5], k_skew=[0.0, 0.0, 0.1])
    assert (back.rc, back.k_normal[0]) == (10.0, 0.0)


def test_exact_cfbend_input_file_geometry(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.write_text(
        "lattice.elements = geo_ecf_b\n"
        "geo_ecf_b.type = cfbend_exact\ngeo_ecf_b.ds = 0.5\n"
        "geo_ecf_b.k_normal = 0.0 0.3\ngeo_ecf_b.k_skew = 0.0 0.0\n"
        "geo_ecf_b.B = 0.5\n"
    )

    sim = ImpactX()
    sim.load_inputs_file(str(inputs))
    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.init_grids()
    sim.init_lattice_elements_from_inputs()

    bend = sim.lattice[0]
    assert (bend.rc, bend.phi, bend.B) == (None, None, 0.5)
    sim.finalize()
