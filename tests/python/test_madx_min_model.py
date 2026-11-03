#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Tests for the ``min_model`` fidelity floor of the MAD-X -> ImpactX translation.

ImpactX implements most elements at up to three levels of fidelity: linear,
paraxial (``Chr*``) and exact (``Exact*``), see ``docs/source/theory/assumptions.rst``.
The MAD-X translator picks the cheapest ImpactX element that represents the
MAD-X input. ``min_model`` raises the *lower gate* of that choice without
otherwise changing it. Three properties are asserted here:

1. A floor never lowers a tier. An element that already needs a richer model,
   such as a thick octupole that only exists as ``ExactMultipole``, keeps it at
   ``min_model="linear"``.
2. Where a tier is not implemented, the floor rounds *up*. ImpactX has no
   ``ChrSbend``, so ``min_model="paraxial"`` on a plain SBEND yields the exact
   model rather than falling back to the linear one.
3. Where no model reaches the floor at all, the translation still succeeds and
   warns once. SOLENOID has no exact model.

References:
    - https://impactx.readthedocs.io/en/latest/usage/python.html#impactx.elements.KnownElementsList.load_file
    - MAD-X manual: https://madx.web.cern.ch/webguide/manual.html
"""

import math
import warnings

import numpy as np
import pytest

from impactx import Config, RefPart, elements
from impactx.element_models import MODEL_TIERS, select_model
from impactx.madx_to_impactx import (
    MADXImpactXTranslatorWarning,
    lattice,
    read_beam,
    read_lattice,
)
from impactx.MADXParser import MADXInputError


def _translate(elems, **kwargs):
    """Translate parsed MAD-X element dicts through lattice(), muting warnings."""
    if isinstance(elems, dict):
        elems = [elems]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return lattice(elems, **kwargs)


def _types(beamline):
    """Element class names of a translated beamline."""
    return [type(e).__name__ for e in beamline]


def _only(beamline, cls):
    """Return the single element of type ``cls`` in a translated beamline."""
    matches = [e for e in beamline if isinstance(e, cls)]
    assert len(matches) == 1, _types(beamline)
    return matches[0]


# ---------------------------------------------------------------------------
# The tier picker itself
# ---------------------------------------------------------------------------


def test_select_model_picks_cheapest_at_or_above_floor():
    """The cheapest implemented tier that satisfies the floor wins."""
    full = {tier: tier for tier in MODEL_TIERS}
    assert select_model(full, "linear") == ("linear", "linear")
    assert select_model(full, "paraxial") == ("paraxial", "paraxial")
    assert select_model(full, "exact") == ("exact", "exact")


def test_select_model_rounds_up_through_a_missing_tier():
    """A family without a paraxial model serves a paraxial floor from its exact one."""
    holed = {"linear": "linear", "exact": "exact"}
    assert select_model(holed, "linear") == ("linear", "linear")
    assert select_model(holed, "paraxial") == ("exact", "exact")
    assert select_model(holed, "exact") == ("exact", "exact")


def test_select_model_falls_back_below_an_unreachable_floor():
    """With nothing at or above the floor, the most faithful model available is used."""
    linear_only = {"linear": "linear"}
    assert select_model(linear_only, "exact") == ("linear", "linear")


# ---------------------------------------------------------------------------
# Element families with a full or partial tier ladder
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "min_model, expected",
    [("linear", "Drift"), ("paraxial", "ChrDrift"), ("exact", "ExactDrift")],
)
def test_drift_tiers(min_model, expected):
    """DRIFT has all three tiers, so the floor is met exactly."""
    beamline = _translate(
        {"name": "d1", "type": "drift", "l": 0.5}, min_model=min_model
    )
    assert _types(beamline) == [expected]
    assert beamline[0].ds == pytest.approx(0.5)


@pytest.mark.parametrize(
    "min_model, expected",
    [("linear", "Quad"), ("paraxial", "ChrQuad"), ("exact", "ExactQuad")],
)
def test_quadrupole_tiers(min_model, expected):
    """QUADRUPOLE has all three tiers. Strength and rotation survive the promotion."""
    beamline = _translate(
        {"name": "q1", "type": "quadrupole", "l": 0.3, "k1": 2.0, "tilt": 0.25},
        min_model=min_model,
    )
    assert _types(beamline) == [expected]
    quad = beamline[0]
    assert quad.ds == pytest.approx(0.3)
    assert quad.k == pytest.approx(2.0)
    assert quad.rotation == pytest.approx(0.25 * 180.0 / math.pi)
    if min_model != "linear":
        # MAD-X k convention: k in m^(-2), which plain Quad always assumes
        assert quad.unit == 0


@pytest.mark.parametrize(
    "min_model, expected",
    [("linear", "Sbend"), ("paraxial", "ExactSbend"), ("exact", "ExactSbend")],
)
def test_sbend_tiers_round_up(min_model, expected):
    """There is no ChrSbend: a paraxial floor rounds up to the exact bend."""
    ds, angle = 2.0, 0.1
    beamline = _translate(
        {"name": "b1", "type": "sbend", "l": ds, "angle": angle},
        min_model=min_model,
    )
    assert _types(beamline) == [expected]
    bend = beamline[0]
    assert bend.ds == pytest.approx(ds)
    if expected == "Sbend":
        assert bend.to_dict()["rc"] == pytest.approx(ds / angle)
    else:
        # note: the ExactSbend constructor takes phi in degrees, but the
        # property returns radians (ImpactX issue #1367)
        assert bend.phi == pytest.approx(angle)


@pytest.mark.parametrize(
    "min_model, expected",
    [("linear", "CFbend"), ("paraxial", "ExactCFbend"), ("exact", "ExactCFbend")],
)
def test_combined_function_bend_tiers_round_up(min_model, expected):
    """A bend with K1 uses the combined-function rung of the ladder."""
    ds, angle, k1 = 2.0, 0.1, 0.5
    beamline = _translate(
        {"name": "b1", "type": "sbend", "l": ds, "angle": angle, "k1": k1},
        min_model=min_model,
    )
    assert _types(beamline) == [expected]
    fields = beamline[0].to_dict()
    if expected == "CFbend":
        assert fields["rc"] == pytest.approx(ds / angle)
        assert fields["k"] == pytest.approx(k1)
    else:
        assert fields["k_normal"][0] == pytest.approx(angle / ds)  # curvature 1/rc
        assert fields["k_normal"][1] == pytest.approx(k1)


@pytest.mark.parametrize(
    "min_model, expected",
    [
        ("linear", "linear_slac75"),
        ("paraxial", "nonlinear_hwang_lee"),
        ("exact", "nonlinear_hwang_lee"),
    ],
)
def test_dipedge_fringe_model_follows_min_model(min_model, expected):
    """DipEdge switches its fringe model by argument. "nonlinear_hwang_lee" is its exact tier."""
    beamline = _translate(
        {
            "name": "de",
            "type": "dipedge",
            "h": 0.0966,
            "e1": 0.0483,
            "hgap": 0.01,
            "fint": 0.5,
        },
        min_model=min_model,
    )
    edge = _only(beamline, elements.DipEdge)
    assert edge.model == expected
    # physical parameters are untouched by the model choice
    assert edge.psi == pytest.approx(0.0483)
    assert edge.rc == pytest.approx(1.0 / 0.0966)


def test_bend_edges_follow_min_model():
    """The DipEdges a bend emits at its faces honor the floor, too."""
    elem = {
        "name": "b1",
        "type": "sbend",
        "l": 2.0,
        "angle": 0.1,
        "e1": 0.05,
        "e2": 0.05,
        "hgap": 0.01,
        "fint": 0.5,
    }
    beamline = _translate(dict(elem), min_model="linear")
    assert _types(beamline) == ["DipEdge", "Sbend", "DipEdge"]
    assert [e.model for e in beamline if isinstance(e, elements.DipEdge)] == [
        "linear_slac75",
        "linear_slac75",
    ]

    beamline = _translate(dict(elem), min_model="exact")
    assert _types(beamline) == ["DipEdge", "ExactSbend", "DipEdge"]
    assert [e.model for e in beamline if isinstance(e, elements.DipEdge)] == [
        "nonlinear_hwang_lee",
        "nonlinear_hwang_lee",
    ]


def test_synthetic_drifts_follow_min_model():
    """Drifts the translator adds itself (here: around a thin MONITOR) follow the floor."""
    beamline = _translate(
        {"name": "m1", "type": "monitor", "l": 0.4}, min_model="exact"
    )
    assert _types(beamline) == ["ExactDrift", "BeamMonitor"]


# ---------------------------------------------------------------------------
# A floor never lowers a tier
# ---------------------------------------------------------------------------


def test_skew_quadrupole_stays_exact_at_linear_floor():
    """A skew quadrupole only exists as ExactMultipole. The floor must not downgrade it."""
    beamline = _translate(
        {"name": "q1", "type": "quadrupole", "l": 0.3, "k1": 2.0, "k1s": 0.4},
        min_model="linear",
    )
    assert _types(beamline) == ["ExactMultipole"]


def test_sextupole_stays_exact_at_linear_floor():
    """ImpactX has no thick linear sextupole. SEXTUPOLE stays exact at any floor."""
    beamline = _translate(
        {"name": "s1", "type": "sextupole", "l": 0.2, "k2": 3.0}, min_model="linear"
    )
    assert _types(beamline) == ["ExactMultipole"]


# ---------------------------------------------------------------------------
# Unreachable floors warn, once
# ---------------------------------------------------------------------------


def test_solenoid_warns_once_when_floor_is_unreachable():
    """Without a reference energy only the linear Sol is available: warn once.

    The paraxial ChrAcc(ez=0) model needs beta*gamma to convert ``ks``, so a
    lattice translated without it keeps ``Sol`` at every floor.
    """
    solenoids = [
        {"name": f"s{i}", "type": "solenoid", "l": 0.5, "ks": 1.0} for i in range(3)
    ]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        beamline = lattice(solenoids, min_model="exact")

    assert _types(beamline) == ["Sol", "Sol", "Sol"]
    floor_warnings = [
        w
        for w in caught
        if issubclass(w.category, MADXImpactXTranslatorWarning)
        and "SOLENOID" in str(w.message)
        and "exact" in str(w.message)
    ]
    assert len(floor_warnings) == 1, [str(w.message) for w in caught]


def test_solenoid_does_not_warn_at_linear_floor():
    """The default floor is reachable everywhere, so it never triggers the warning."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        beamline = lattice([{"name": "s1", "type": "solenoid", "l": 0.5, "ks": 1.0}])

    assert _types(beamline) == ["Sol"]
    assert not [w for w in caught if "no 'linear'" in str(w.message)]


# ---------------------------------------------------------------------------
# Solenoid promotion to the paraxial ChrAcc(ez=0) model
# ---------------------------------------------------------------------------

# beta*gamma of a 2 GeV kinetic-energy proton, the energy used below.
_PROTON_MASS_MeV = 938.27208816
_SOL_BETA_GAMMA = math.sqrt(((2000.0 + _PROTON_MASS_MeV) / _PROTON_MASS_MeV) ** 2 - 1.0)


@pytest.mark.parametrize(
    "min_model,expected",
    [("linear", "Sol"), ("paraxial", "ChrAcc"), ("exact", "ChrAcc")],
)
def test_solenoid_promotes_to_chracc_with_reference_energy(min_model, expected):
    """With beta*gamma known, the paraxial tier is ChrAcc(ez=0); exact rounds down."""
    beamline = _translate(
        {"name": "s1", "type": "solenoid", "l": 0.5, "ks": 1.0},
        min_model=min_model,
        ref_beta_gamma=_SOL_BETA_GAMMA,
    )
    assert _types(beamline) == [expected]


def test_solenoid_chracc_field_follows_rigidity_conversion():
    """Sol.ks is per unit rigidity, ChrAcc.bz is not: bz = ks * beta_gamma.

    Both maps rotate by the same angle for the reference particle, Sol by
    ks/2 * ds and ChrAcc by bz/2 * ds / beta_gamma.
    """
    ks = 1.7
    beamline = _translate(
        {"name": "s1", "type": "solenoid", "l": 0.5, "ks": ks},
        min_model="paraxial",
        ref_beta_gamma=_SOL_BETA_GAMMA,
    )
    sol = _only(beamline, elements.ChrAcc)
    assert sol.ez == 0.0
    assert sol.bz == pytest.approx(ks * _SOL_BETA_GAMMA)


@pytest.mark.parametrize("ks", [1.7, -0.9])
def test_solenoid_promotion_preserves_the_reference_map(ks):
    """The promoted ChrAcc is the same transport map as the Sol it replaces.

    ``ChrAcc`` adds the chromatic dependence on top, so the two agree exactly
    only for the reference particle, which is what the 6x6 map describes.

    ``ks=0`` is left out: ``Sol`` divides by its own strength when building the
    map and returns NaN there, independent of this translation.
    """
    ref = RefPart()
    ref.set_species("proton").set_kin_energy_MeV(2000.0)

    common = {"name": "s1", "type": "solenoid", "l": 0.5, "ks": ks}
    sol = _only(_translate(common, min_model="linear"), elements.Sol)
    acc = _only(
        _translate(common, min_model="paraxial", ref_beta_gamma=ref.beta_gamma),
        elements.ChrAcc,
    )

    linear_map = np.array(sol.transfer_map(ref)).reshape(6, 6)
    paraxial_map = np.array(acc.transfer_map(ref)).reshape(6, 6)
    atol = 1.0e-14 if Config.precision != "SINGLE" else 1.0e-6
    assert np.allclose(linear_map, paraxial_map, rtol=0.0, atol=atol)


def test_solenoid_promotion_warns_about_energy_pinning():
    """ChrAcc.bz is fixed at translation time, so the reader says so, once."""
    solenoids = [
        {"name": f"s{i}", "type": "solenoid", "l": 0.5, "ks": 1.0} for i in range(3)
    ]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        beamline = lattice(
            solenoids, min_model="paraxial", ref_beta_gamma=_SOL_BETA_GAMMA
        )

    assert _types(beamline) == ["ChrAcc"] * 3
    pinning = [w for w in caught if "pinned to the MAD-X BEAM energy" in str(w.message)]
    assert len(pinning) == 1, [str(w.message) for w in caught]
    # The paraxial floor is reachable now, so no "no model available" warning.
    assert not [w for w in caught if "model available" in str(w.message)]


def test_solenoid_exact_floor_warns_but_still_promotes():
    """There is no ExactAcc yet, so an exact floor rounds down to ChrAcc and warns."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        beamline = lattice(
            [{"name": "s1", "type": "solenoid", "l": 0.5, "ks": 1.0}],
            min_model="exact",
            ref_beta_gamma=_SOL_BETA_GAMMA,
        )

    assert _types(beamline) == ["ChrAcc"]
    floor = [
        w
        for w in caught
        if issubclass(w.category, MADXImpactXTranslatorWarning)
        and "SOLENOID has no 'exact'" in str(w.message)
        and "'paraxial' model instead" in str(w.message)
    ]
    assert len(floor) == 1, [str(w.message) for w in caught]


@pytest.mark.parametrize(
    "beam_energy",
    # 2 GeV kinetic energy as ENERGY and as PC
    ["energy=2.93827208816", "pc=sqrt(2 * (2 + 2 * pmass))"],
)
def test_solenoid_promotion_end_to_end_from_madx_file(tmp_path, beam_energy):
    """read_lattice derives beta*gamma from BEAM and promotes the solenoid."""
    madx_file = tmp_path / "sol.madx"
    madx_file.write_text(
        f"beam, particle=proton, {beam_energy};\n"
        "s1: solenoid, l=0.5, ks=1.0;\n"
        "cell: line=(s1);\n"
        "use, sequence=cell;\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        linear = read_lattice(str(madx_file), min_model="linear")
        paraxial = read_lattice(str(madx_file), min_model="paraxial")

    assert _types(linear) == ["Sol"]
    assert _types(paraxial) == ["ChrAcc"]
    assert _only(paraxial, elements.ChrAcc).bz == pytest.approx(
        1.0 * _SOL_BETA_GAMMA, rel=1e-6
    )


@pytest.mark.parametrize(
    "beam_commands, beta_gamma",
    [
        ("beam, particle=ion, mass=200, charge=82, pc=300;", 1.5),
        ("beam, particle=ion, mass=200, charge=82, gamma=5;", math.sqrt(24.0)),
        ("beam, particle=proton, mass=2, pc=3;", 3.0 / 0.93827208816),
        (
            "beam, mass=2, charge=1, pc=3;\nbeam, particle=proton, pc=2;",
            2.0 / 0.93827208816,
        ),
        # a stated momentum is kept, also below the resolution of the total energy
        ("beam, particle=proton, pc=1e-9;", 1.0e-9 / 0.93827208816),
    ],
)
def test_solenoid_promotion_uses_the_particle_of_the_beam(
    tmp_path, beam_commands, beta_gamma
):
    """beta*gamma uses the same mass as the energy, which MAD-X resolves per BEAM."""
    madx_file = tmp_path / "sol.madx"
    madx_file.write_text(
        f"{beam_commands}\n"
        "s1: solenoid, l=0.5, ks=1.0;\n"
        "cell: line=(s1);\n"
        "use, sequence=cell;\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        paraxial = read_lattice(str(madx_file), min_model="paraxial")

    rel = 1.0e-14 if Config.precision != "SINGLE" else 1.0e-6
    assert _only(paraxial, elements.ChrAcc).bz == pytest.approx(beta_gamma, rel=rel)


@pytest.mark.parametrize(
    "beam_command, charge_qe, mass_MeV, beta_gamma",
    [
        ("beam;", 1.0, 0.51099895, 0.99999986944002806 / 0.51099895e-3),
        ("beam, particle=ion, mass=200, charge=82, pc=300;", 82.0, 200.0e3, 1.5),
        (
            "beam, particle=ion, mass=200, charge=82, gamma=5;",
            82.0,
            200.0e3,
            math.sqrt(24.0),
        ),
        ("beam, particle=proton, pc=2;", 1.0, 938.27208816, 2.0 / 0.93827208816),
        ("beam, particle=proton, beta=0;", 1.0, 938.27208816, 0.0),
        ("beam, particle=proton, brho=0;", 1.0, 938.27208816, 0.0),
        ("beam, particle=posmuon, gamma=5;", 1.0, 105.6583755, math.sqrt(24.0)),
        ("beam, particle=antiproton, gamma=5;", -1.0, 938.27208816, math.sqrt(24.0)),
    ],
)
def test_read_beam(tmp_path, beam_command, charge_qe, mass_MeV, beta_gamma):
    """RefPart.load_file takes the particle and energy MAD-X resolves, also at rest."""
    madx_file = tmp_path / "beam.madx"
    madx_file.write_text(f"{beam_command}\n")
    ref = RefPart()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        read_beam(ref, str(madx_file))

    is_double = Config.precision != "SINGLE"
    assert ref.charge_qe == pytest.approx(charge_qe, rel=1e-14 if is_double else 1e-6)
    assert ref.mass_MeV == pytest.approx(mass_MeV, rel=1e-14 if is_double else 1e-6)
    assert ref.beta_gamma == pytest.approx(
        beta_gamma, rel=1e-9 if is_double else 1e-6, abs=1e-12
    )


def test_read_beam_low_energy(tmp_path):
    """RefPart.load_file keeps a low kinetic energy."""
    madx_file = tmp_path / "beam.madx"
    madx_file.write_text("beam, particle=proton, beta=1e-3;\n")
    ref = RefPart()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        read_beam(ref, str(madx_file))

    # single precision: the reference `pt` near -1 cannot resolve a beta*gamma of 1e-3
    rel = 1e-9 if Config.precision != "SINGLE" else 3e-2
    assert ref.beta_gamma == pytest.approx(1.0e-3 / math.sqrt(1.0 - 1.0e-6), rel=rel)


@pytest.mark.parametrize("beam_command", ["", "beam, mass=-1, charge=1;"])
def test_rf_voltage_needs_the_reference_mass(tmp_path, beam_command):
    """A voltage cannot be normalized without a mass, unlike a cavity that is off."""
    madx_file = tmp_path / "rf.madx"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        madx_file.write_text(
            f"{beam_command}\nrf: rfcavity, volt=5, freq=100;\ns1: line=(rf);\n"
        )
        with pytest.raises(MADXInputError, match="RFCAVITY with VOLT=5"):
            read_lattice(str(madx_file))

        madx_file.write_text(
            f"{beam_command}\nrf: rfcavity, volt=0, freq=100;\ns1: line=(rf);\n"
        )
        assert _only(read_lattice(str(madx_file)), elements.ShortRF).V == 0.0


@pytest.mark.parametrize(
    "definition",
    [
        "drift, l=sqrt(-1)",
        "drift, l=(-1)^0.5",
        "quadrupole, l=1, k1=sqrt(-1)",
        "sbend, l=1, angle=sqrt(-1)",
        "drift, l=exp(1000)",
    ],
)
def test_element_in_use_with_a_non_finite_attribute_raises(tmp_path, definition):
    """The translation would treat NaN like zero, and drop or replace the element."""
    madx_file = tmp_path / "element.madx"
    madx_file.write_text(
        f"beam, particle=proton, energy=2;\ne1: {definition};\ns1: line=(e1);\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(MADXInputError, match="'e1'.*not finite"):
            read_lattice(str(madx_file))


def test_element_not_in_use_with_a_non_finite_attribute(tmp_path):
    """An element that is only defined does not matter, as in MAD-X."""
    madx_file = tmp_path / "element.madx"
    madx_file.write_text(
        "beam, particle=proton, energy=2;\n"
        "unused: drift, l=sqrt(-1);\n"
        "d1: drift, l=1;\n"
        "s1: line=(d1);\n"
        "use, sequence=s1;\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert _types(read_lattice(str(madx_file))) == ["Drift"]


@pytest.mark.parametrize("argument", ["line", "sequence"])
def test_read_lattice_uses_the_beam_of_the_requested_line(tmp_path, argument):
    """The requested line counts, not the one of USE, which has no beam here."""
    madx_file = tmp_path / "lines.madx"
    madx_file.write_text(
        "d1: drift, l=1;\n"
        "a: line=(d1);\n"
        "b: line=(d1);\n"
        "beam, sequence=b, particle=proton, energy=2;\n"
        "use, period=a;\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert _types(read_lattice(str(madx_file), **{argument: "b"})) == ["Drift"]


@pytest.mark.parametrize(
    "use_commands, name",
    [
        ("use, sequence=food;", "food"),
        ("use, sequence=fodo;\nuse, sequence=food;", "food"),
        ('use, sequence="FODO";', "FODO"),
    ],
)
def test_use_of_an_unknown_sequence_raises_in_the_lattice_import(
    tmp_path, use_commands, name
):
    """Stricter than MAD-X, which warns: a typo in USE must not load another line."""
    madx_file = tmp_path / "use.madx"
    madx_file.write_text(
        "d1: drift, l=1;\n"
        "fodo: line=(d1);\n"
        "beam, particle=proton, energy=2;\n"
        f"{use_commands}\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(MADXInputError, match=f"sequence '{name}' is not defined"):
            read_lattice(str(madx_file))

        # a requested line does not depend on it, and neither does the reference particle
        assert _types(read_lattice(str(madx_file), line="fodo")) == ["Drift"]
        ref = RefPart()
        read_beam(ref, str(madx_file))
    rel = 1.0e-14 if Config.precision != "SINGLE" else 1.0e-6
    assert ref.mass_MeV == pytest.approx(938.27208816, rel=rel)


def test_use_of_a_known_sequence_replaces_an_unknown_one(tmp_path):
    """Only the last USE counts, also for a sequence that an earlier one could not find."""
    madx_file = tmp_path / "use.madx"
    madx_file.write_text(
        "d1: drift, l=1;\n"
        "fodo: line=(d1);\n"
        "beam, particle=proton, energy=2;\n"
        "use, sequence=food;\n"
        "use, sequence=fodo;\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert _types(read_lattice(str(madx_file))) == ["Drift"]


def test_read_lattice_uses_the_beam_of_the_inferred_sequence(tmp_path):
    """Without USE, the lattice import selects the line, and that line's beam."""
    madx_file = tmp_path / "sol.madx"
    madx_file.write_text(
        "beam, particle=electron, energy=2;\n"
        "beam, sequence=cell, particle=proton, pc=2;\n"
        "s1: solenoid, l=0.5, ks=1.0;\n"
        "cell: line=(s1);\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        paraxial = read_lattice(str(madx_file), min_model="paraxial")

    rel = 1.0e-14 if Config.precision != "SINGLE" else 1.0e-6
    assert _only(paraxial, elements.ChrAcc).bz == pytest.approx(
        2.0 / 0.93827208816, rel=rel
    )


def test_read_lattice_without_a_beam_command(tmp_path):
    """MAD-X stops at USE without BEAM, but reads a file of definitions."""
    definitions = "d1: drift, l=0.5;\ncell: line=(d1);\n"
    madx_file = tmp_path / "no_beam.madx"

    madx_file.write_text(definitions)
    with pytest.warns(MADXImpactXTranslatorWarning, match="no BEAM command"):
        assert _types(read_lattice(str(madx_file))) == ["Drift"]

    madx_file.write_text(f"{definitions}use, sequence=cell;\n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(MADXInputError, match="no BEAM command"):
            read_lattice(str(madx_file))

    # an empty BEAM is enough
    madx_file.write_text(f"beam;\n{definitions}use, sequence=cell;\n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert _types(read_lattice(str(madx_file))) == ["Drift"]


def test_read_beam_without_a_beam_command(tmp_path):
    """There is no reference particle to load from a file without BEAM."""
    madx_file = tmp_path / "no_beam.madx"
    madx_file.write_text("d1: drift, l=0.5;\ncell: line=(d1);\n")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(MADXInputError, match="no BEAM command"):
            read_beam(RefPart(), str(madx_file))


@pytest.mark.parametrize(
    "beam_command",
    [
        "beam, mass=-1, charge=1;",
        "beam, mass=-1, charge=1, energy=0.5;",
        "beam, particle=proton, beta=-1;",
        "beam, particle=proton, beta=-2;",
        "beam, particle=proton, pc=sqrt(-1);",
        "beam, particle=proton, gamma=asin(2);",
        "beam, mass=sqrt(-1), charge=1, energy=3;",
        "beam, mass=2, charge=sqrt(-1), energy=3;",
    ],
)
def test_beam_that_is_no_reference_particle(tmp_path, beam_command):
    """MAD-X reads these, and so does the lattice import: only RefPart cannot use them."""
    madx_file = tmp_path / "beam.madx"
    madx_file.write_text(
        f"{beam_command}\nd1: drift, l=0.5;\ncell: line=(d1);\nuse, sequence=cell;\n"
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(MADXInputError, match="no valid reference particle"):
            read_beam(RefPart(), str(madx_file))

    with pytest.warns(
        MADXImpactXTranslatorWarning, match="no valid reference particle"
    ):
        assert _types(read_lattice(str(madx_file))) == ["Drift"]


def test_solenoid_not_promoted_without_a_stated_beam_energy(tmp_path):
    """MAD-X defaults ENERGY to 1 GeV: do not translate a rigidity against it.

    Promoting here would bake a wrong field strength into the lattice, so the
    reader keeps the momentum-independent Sol and says why.
    """
    madx_file = tmp_path / "sol_no_energy.madx"
    madx_file.write_text(
        "beam, particle=proton;\n"  # no ENERGY: MAD-X assumes 1 GeV
        "s1: solenoid, l=0.5, ks=1.0;\n"
        "cell: line=(s1);\n"
        "use, sequence=cell;\n"
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        beamline = read_lattice(str(madx_file), min_model="paraxial")

    assert _types(beamline) == ["Sol"]
    explained = [w for w in caught if "does not provide" in str(w.message)]
    assert len(explained) == 1, [str(w.message) for w in caught]


# ---------------------------------------------------------------------------
# Validation and end-to-end plumbing
# ---------------------------------------------------------------------------


def test_invalid_min_model_raises():
    with pytest.raises(ValueError, match="min_model"):
        lattice([{"name": "d1", "type": "drift", "l": 0.5}], min_model="quadratic")


def test_read_lattice_forwards_min_model(tmp_path):
    """min_model reaches the translator through the file-reading entry point."""
    deck = tmp_path / "fodo.madx"
    deck.write_text(
        "BEAM, PARTICLE=ELECTRON, ENERGY=5.0;\n"
        "D1: DRIFT, L=0.25;\n"
        "Q1: QUADRUPOLE, L=0.5, K1=1.0;\n"
        "MYLINE: LINE=(D1, Q1, D1);\n"
        "USE, SEQUENCE=MYLINE;\n"
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        default = read_lattice(str(deck))
        exact = read_lattice(str(deck), min_model="exact")

    assert _types(default) == ["Drift", "Quad", "Drift"]
    assert _types(exact) == ["ExactDrift", "ExactQuad", "ExactDrift"]


def test_load_file_forwards_min_model(tmp_path):
    """The public KnownElementsList.load_file entry point forwards min_model."""
    deck = tmp_path / "fodo.madx"
    deck.write_text(
        "BEAM, PARTICLE=ELECTRON, ENERGY=5.0;\n"
        "D1: DRIFT, L=0.25;\n"
        "Q1: QUADRUPOLE, L=0.5, K1=1.0;\n"
        "MYLINE: LINE=(D1, Q1, D1);\n"
        "USE, SEQUENCE=MYLINE;\n"
    )

    lat = elements.KnownElementsList()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lat.load_file(str(deck), nslice=2, min_model="paraxial")

    assert _types(lat) == ["ChrDrift", "ChrQuad", "ChrDrift"]
    assert lat[1].nslice == 2


def test_load_file_rejects_invalid_min_model(tmp_path):
    deck = tmp_path / "fodo.madx"
    deck.write_text(
        "BEAM, PARTICLE=ELECTRON, ENERGY=5.0;\n"
        "D1: DRIFT, L=0.25;\n"
        "MYLINE: LINE=(D1);\n"
        "USE, SEQUENCE=MYLINE;\n"
    )

    lat = elements.KnownElementsList()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="min_model"):
            lat.load_file(str(deck), min_model="chromatic")


# ---------------------------------------------------------------------------
# PALS reader (same floor, same vocabulary)
# ---------------------------------------------------------------------------


def test_pals_min_model():
    """The PALS reader shares the tier tables. Its default output is unchanged."""
    pals = pytest.importorskip("pals")

    beamline = pals.BeamLine(
        name="fodo",
        line=[
            pals.Drift(name="d1", length=0.25),
            pals.Quadrupole(
                name="q1",
                length=0.5,
                MagneticMultipoleP=pals.MagneticMultipoleParameters(Kn1=1.0),
            ),
        ],
    )

    from impactx.pals_to_impactx import read_lattice as pals_read_lattice

    # the linear `Quad` has no `unit` argument yet (ImpactX issue #798), so the
    # cheapest quadrupole model available to PALS is the paraxial one
    assert _types(pals_read_lattice(beamline)) == ["Drift", "ChrQuad"]
    assert _types(pals_read_lattice(beamline, min_model="exact")) == [
        "ExactDrift",
        "ExactQuad",
    ]
