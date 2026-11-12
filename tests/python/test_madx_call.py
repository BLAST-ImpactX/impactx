#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Tests for the MAD-X parser: CALL command, line continuation, and general syntax.

References:
    - https://github.com/BLAST-ImpactX/impactx/issues/512
    - http://madx.web.cern.ch/madx/webguide/manual.html (Chapter 5.2)
    - MAD-X source: src/mad_write.c, src/mad_cmd.c, src/mad_parse.c
"""

import math
import re
import warnings

import pytest

from impactx.MADXParser import MADXInputError, MADXInputWarning, MADXParser


@pytest.fixture
def tmp_madx(tmp_path):
    """Helper to write temporary MAD-X files and return their paths."""

    def _write(filename, content):
        filepath = tmp_path / filename
        filepath.parent.mkdir(parents=True, exist_ok=True)
        filepath.write_text(content)
        return str(filepath)

    return _write


def test_call_file(tmp_madx):
    """Test basic CALL, FILE='...' command (MAD-X syntax)."""
    elements_file = tmp_madx(
        "elements.madx",
        """\
D1: DRIFT, L=0.25;
QF: QUADRUPOLE, L=1.0, K1=1.0;
QD: QUADRUPOLE, L=1.0, K1=-1.0;
""",
    )

    main_file = tmp_madx(
        "main.madx",
        f"""\
BEAM, PARTICLE=ELECTRON, ENERGY=2.0;
CALL, FILE="{elements_file}";
FODO: LINE=(D1, QF, D1, QD);
USE, SEQUENCE=FODO;
""",
    )

    parser = MADXParser()
    parser.parse(main_file)

    beamline = parser.getBeamline()
    assert len(beamline) == 4
    assert beamline[0]["type"] == "drift"
    assert beamline[0]["l"] == 0.25
    assert beamline[1]["type"] == "quadrupole"
    assert beamline[1]["k1"] == 1.0
    assert beamline[2]["type"] == "drift"
    assert beamline[3]["type"] == "quadrupole"
    assert beamline[3]["k1"] == -1.0
    assert parser.getParticle() == "electron"


def test_call_filename(tmp_madx):
    """Test CALL, FILENAME='...' command (MAD-8 syntax)."""
    elements_file = tmp_madx(
        "elements.madx",
        """\
D1: DRIFT, L=0.5;
""",
    )

    main_file = tmp_madx(
        "main.madx",
        f"""\
BEAM, PARTICLE=PROTON, ENERGY=7000.0;
CALL, FILENAME="{elements_file}";
MYLINE: LINE=(D1);
USE, SEQUENCE=MYLINE;
""",
    )

    parser = MADXParser()
    parser.parse(main_file)

    beamline = parser.getBeamline()
    assert len(beamline) == 1
    assert beamline[0]["type"] == "drift"
    assert beamline[0]["l"] == 0.5
    assert parser.getParticle() == "proton"


def test_call_nested(tmp_madx):
    """Test nested CALL statements (CALL within a CALL'd file)."""
    inner_file = tmp_madx(
        "inner.madx",
        """\
QF: QUADRUPOLE, L=1.0, K1=2.5;
""",
    )

    outer_file = tmp_madx(
        "outer.madx",
        f"""\
D1: DRIFT, L=0.3;
CALL, FILE="{inner_file}";
""",
    )

    main_file = tmp_madx(
        "main.madx",
        f"""\
BEAM, PARTICLE=ELECTRON, ENERGY=3.0;
CALL, FILE="{outer_file}";
MYLINE: LINE=(D1, QF, D1);
USE, SEQUENCE=MYLINE;
""",
    )

    parser = MADXParser()
    parser.parse(main_file)

    beamline = parser.getBeamline()
    assert len(beamline) == 3
    assert beamline[0]["type"] == "drift"
    assert beamline[0]["l"] == 0.3
    assert beamline[1]["type"] == "quadrupole"
    assert beamline[1]["k1"] == 2.5
    assert beamline[2]["type"] == "drift"


def test_call_with_return(tmp_madx):
    """Test RETURN statement stops reading the called file."""
    called_file = tmp_madx(
        "partial.madx",
        """\
D1: DRIFT, L=0.1;
RETURN;
! Everything below should be ignored
QF: QUADRUPOLE, L=1.0, K1=99.0;
""",
    )

    main_file = tmp_madx(
        "main.madx",
        f"""\
BEAM, PARTICLE=ELECTRON, ENERGY=1.0;
CALL, FILE="{called_file}";
MYLINE: LINE=(D1);
USE, SEQUENCE=MYLINE;
""",
    )

    parser = MADXParser()
    parser.parse(main_file)

    beamline = parser.getBeamline()
    assert len(beamline) == 1
    assert beamline[0]["type"] == "drift"
    # QF should NOT have been parsed
    assert parser.getElement("qf") is None


def test_call_variables_shared(tmp_madx):
    """Test that variables defined in called files are visible in the main file."""
    vars_file = tmp_madx(
        "variables.madx",
        """\
mylen = 0.75;
myk1 = 3.14;
""",
    )

    main_file = tmp_madx(
        "main.madx",
        f"""\
BEAM, PARTICLE=ELECTRON, ENERGY=1.0;
CALL, FILE="{vars_file}";
D1: DRIFT, L=mylen;
QF: QUADRUPOLE, L=mylen, K1=myk1;
MYLINE: LINE=(D1, QF);
USE, SEQUENCE=MYLINE;
""",
    )

    parser = MADXParser()
    parser.parse(main_file)

    beamline = parser.getBeamline()
    assert len(beamline) == 2
    assert beamline[0]["l"] == 0.75
    assert beamline[1]["l"] == 0.75
    assert beamline[1]["k1"] == 3.14


def test_call_relative_path(tmp_madx):
    """Test CALL with a file referenced relative to the calling file."""
    # Create a subdirectory structure
    _ = tmp_madx(
        "subdir/elements.madx",
        """\
D1: DRIFT, L=0.42;
""",
    )

    # Reference relative to the main file's directory
    main_file = tmp_madx(
        "main.madx",
        """\
BEAM, PARTICLE=ELECTRON, ENERGY=1.0;
CALL, FILE="subdir/elements.madx";
MYLINE: LINE=(D1);
USE, SEQUENCE=MYLINE;
""",
    )

    parser = MADXParser()
    parser.parse(main_file)

    beamline = parser.getBeamline()
    assert len(beamline) == 1
    assert beamline[0]["l"] == 0.42


def test_call_file_not_found(tmp_madx):
    """Test that CALL with a nonexistent file raises an error."""
    main_file = tmp_madx(
        "main.madx",
        """\
CALL, FILE="does_not_exist.madx";
""",
    )

    parser = MADXParser()
    with pytest.raises((FileNotFoundError, MADXInputError)):
        parser.parse(main_file)


def test_call_missing_file_attribute(tmp_madx):
    """Test that CALL without FILE= raises an error."""
    main_file = tmp_madx(
        "main.madx",
        """\
CALL;
""",
    )

    parser = MADXParser()
    with pytest.raises(MADXInputError):
        parser.parse(main_file)


def test_call_multiple_files(tmp_madx):
    """Test multiple CALL commands in sequence."""
    file_a = tmp_madx(
        "beam.madx",
        """\
BEAM, PARTICLE=ELECTRON, ENERGY=5.0;
""",
    )

    file_b = tmp_madx(
        "elements.madx",
        """\
D1: DRIFT, L=1.0;
QF: QUADRUPOLE, L=0.5, K1=1.2;
""",
    )

    file_c = tmp_madx(
        "lattice.madx",
        """\
MYLINE: LINE=(D1, QF, D1);
USE, SEQUENCE=MYLINE;
""",
    )

    main_file = tmp_madx(
        "main.madx",
        f"""\
CALL, FILE="{file_a}";
CALL, FILE="{file_b}";
CALL, FILE="{file_c}";
""",
    )

    parser = MADXParser()
    parser.parse(main_file)

    assert parser.getParticle() == "electron"
    assert parser.getEtot() == 5.0

    beamline = parser.getBeamline()
    assert len(beamline) == 3
    assert beamline[0]["type"] == "drift"
    assert beamline[1]["type"] == "quadrupole"
    assert beamline[2]["type"] == "drift"


# ---- Tests for new MAD-X parser features (2026-03) ----


def test_while_loop():
    """Test WHILE loop with a counter variable."""
    parser = MADXParser()
    parser.parse_string(
        """\
n = 0;
WHILE (n < 5) {
    n = n + 1;
};
"""
    )
    assert parser.context.get_variable("n") == 5.0


def test_if_elseif_logical_operators():
    """Test IF/ELSEIF with &&, ||, and <> operators."""
    parser = MADXParser()
    parser.parse_string(
        """\
a = 3;
b = 5;
result = 0;
IF (a > 1 && b < 3) {
    result = 1;
} ELSEIF (a < 1 || b > 10) {
    result = 2;
} ELSEIF (a <> 3) {
    result = 3;
} ELSE {
    result = 4;
};
"""
    )
    # a>1 && b<3 is false; a<1||b>10 is false; a<>3 is false => ELSE
    assert parser.context.get_variable("result") == 4.0


def test_exec_macro():
    """Test MACRO definition and EXEC with textual substitution.

    In MAD-X, the $ prefix appears on EXEC arguments, not in the macro body.
    The param name in the body is replaced by the argument value.
    """
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=7000.0;
D1: DRIFT, L=1.0;
QF1: QUADRUPOLE, L=0.5, K1=1.2;
QD1: QUADRUPOLE, L=0.5, K1=-1.2;
MakeFODO(QFName, QDName): MACRO = {
    CELL: LINE=(QFName, D1, QDName, D1);
};
EXEC, MakeFODO($QF1, $QD1);
USE, SEQUENCE=CELL;
"""
    )
    beamline = parser.getBeamline()
    assert len(beamline) == 4
    assert beamline[0]["type"] == "quadrupole"
    assert beamline[0]["k1"] == 1.2
    assert beamline[1]["type"] == "drift"
    assert beamline[2]["type"] == "quadrupole"
    assert beamline[2]["k1"] == -1.2
    assert beamline[3]["type"] == "drift"


def test_sequence_refer_and_drifts():
    """Test SEQUENCE with REFER=ENTRY and automatic drift insertion."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=7000.0;
QF: QUADRUPOLE, L=0.5, K1=1.0;
QD: QUADRUPOLE, L=0.5, K1=-1.0;
RING: SEQUENCE, L=4.0, REFER=ENTRY;
    QF, AT=0.0;
    QD, AT=2.0;
ENDSEQUENCE;
USE, SEQUENCE=RING;
"""
    )
    beamline = parser.getBeamline()
    types = [e["type"] for e in beamline]
    # With REFER=ENTRY: QF entry at 0.0 (exit 0.5), gap 1.5, QD entry at 2.0 (exit 2.5)
    # Expect: QF, drift(1.5m), QD
    assert types == ["quadrupole", "drift", "quadrupole"]
    assert beamline[0]["k1"] == 1.0
    assert beamline[2]["k1"] == -1.0
    assert abs(beamline[1]["l"] - 1.5) < 1e-12


def test_multiple_beam_commands():
    """Test multiple BEAM commands keyed by sequence name."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=100.0;
BEAM, PARTICLE=ELECTRON, ENERGY=5.0, SEQUENCE=ERING;
D1: DRIFT, L=1.0;
ERING: LINE=(D1);
USE, SEQUENCE=ERING;
"""
    )
    # After USE ERING, the beam for ERING should be selected
    assert parser.getParticle() == "electron"
    assert parser.getEtot() == 5.0


@pytest.mark.parametrize(
    "particle, mass_name, pc", [("proton", "pmass", 2.0), ("electron", "emass", 1e-3)]
)
def test_beam_energy_parameters(particle, mass_name, pc):
    """PC, GAMMA, BETA and BRHO of the same beam give its total energy."""
    mass = MADXParser().getVariable(mass_name)
    energy = math.hypot(pc, mass)
    for param in (
        f"PC={pc!r}",
        f"GAMMA={energy / mass!r}",
        f"BETA={pc / energy!r}",
        f"BRHO={pc / 0.299792458!r}",
    ):
        parser = MADXParser()
        parser.parse_string(f"BEAM, PARTICLE={particle}, {param};\n")
        assert parser.getEtot() == pytest.approx(energy, rel=1e-14), param


def test_beam_energy_precedence():
    """Within one BEAM command, ENERGY > PC > GAMMA > BETA > BRHO."""
    pmass = 0.93827208816
    for params, energy in (
        ("ENERGY=3, PC=2", 3.0),
        ("BRHO=3, GAMMA=5, PC=2", math.hypot(2.0, pmass)),
        ("BRHO=3, BETA=0.6", 1.25 * pmass),
        ("ENERGY=3, GAMMA=sqrt(-1)", 3.0),  # an ignored value is not evaluated
        ("ENERGY=3, PC=sqrt(-1)", 3.0),
    ):
        parser = MADXParser()
        with pytest.warns(MADXInputWarning, match="ignoring"):
            parser.parse_string(f"BEAM, PARTICLE=PROTON, {params};\n")
        assert parser.getEtot() == pytest.approx(energy, rel=1e-14), params


def test_beam_energy_per_command():
    """Each BEAM command fixes the energy with the current particle mass."""
    proton_pc2 = math.hypot(2.0, 0.93827208816)

    parser = MADXParser()
    parser.parse_string("BEAM, PARTICLE=PROTON, ENERGY=3;\nBEAM, PC=2;\n")
    assert parser.getEtot() == pytest.approx(proton_pc2, rel=1e-14)

    # a new particle keeps the energy
    parser = MADXParser()
    parser.parse_string("BEAM, PARTICLE=PROTON, PC=2;\nBEAM, PARTICLE=ELECTRON;\n")
    assert parser.getParticle() == "electron"
    assert parser.getEtot() == pytest.approx(proton_pc2, rel=1e-14)

    # ... which must exceed its mass
    with pytest.raises(MADXInputError, match="mass"):
        MADXParser().parse_string(
            "BEAM, PARTICLE=ELECTRON, PC=0.01;\nBEAM, PARTICLE=PROTON;\n"
        )
    parser = MADXParser()
    parser.parse_string("BEAM, PARTICLE=PROTON, BETA=0;\n")  # at rest: E = m
    assert parser.getEtot() == 0.93827208816
    with pytest.raises(MADXInputError, match="mass"):
        MADXParser().parse_string("BEAM, PARTICLE=PROTON, BETA=0;\nBEAM, BV=-1;\n")


@pytest.mark.parametrize("param", ["PC=0", "GAMMA=1", "BETA=1", "BRHO=-1"])
def test_beam_energy_invalid_raises(param):
    """Invalid energy parameters raise."""
    with pytest.raises(MADXInputError, match="BEAM"):
        MADXParser().parse_string(f"BEAM, PARTICLE=PROTON, {param};\n")


@pytest.mark.parametrize(
    "script, sequence, energy, mass, charge",
    [
        # ions take MASS and CHARGE, and otherwise keep the charge of the beam
        (
            "BEAM, PARTICLE=ION, MASS=200, CHARGE=82, GAMMA=5;",
            None,
            1000.0,
            200.0,
            82.0,
        ),
        (
            "BEAM, PARTICLE=ION, MASS=200, CHARGE=82, BETA=0.5;",
            None,
            230.94010767585033,
            200.0,
            82.0,
        ),
        (
            "BEAM, PARTICLE=ION, MASS=200, CHARGE=82, BRHO=10;",
            None,
            316.9105523933585,
            200.0,
            82.0,
        ),
        (
            "BEAM, PARTICLE=ION, MASS=200, CHARGE=82, PC=300;",
            None,
            360.5551275463989,
            200.0,
            82.0,
        ),
        ("BEAM, PARTICLE=ION, ENERGY=3;", None, 3.0, 0.93956542052, 1.0),
        (
            "BEAM, PARTICLE=ELECTRON, ENERGY=2; BEAM, PARTICLE=ION, ENERGY=3;",
            None,
            3.0,
            0.93956542052,
            -1.0,
        ),
        (
            "BEAM, PARTICLE=ION, MASS=200, CHARGE=82, ENERGY=300; BEAM, CHARGE=5, BRHO=10;",
            None,
            200.56093536549983,
            200.0,
            5.0,
        ),
        # a known particle replaces any MASS and CHARGE
        (
            "BEAM, MASS=2, CHARGE=1, PC=3; BEAM, PARTICLE=PROTON, PC=2;",
            None,
            2.209152441869987,
            0.93827208816,
            1.0,
        ),
        (
            "BEAM, PARTICLE=PROTON, MASS=2, CHARGE=3, PC=3;",
            None,
            3.1433031211482176,
            0.93827208816,
            1.0,
        ),
        # MASS without PARTICLE replaces the particle; CHARGE alone is not read
        (
            "BEAM, PARTICLE=PROTON, PC=2; BEAM, MASS=2, CHARGE=1, PC=3;",
            None,
            3.605551275463989,
            2.0,
            1.0,
        ),
        (
            "BEAM, MASS=2, CHARGE=1, ENERGY=3; BEAM, CHARGE=5, BRHO=1;",
            None,
            2.0223440651564912,
            2.0,
            1.0,
        ),
        # a name can continue with signs and numbers, and the charge defaults to +1
        (
            "BEAM, PARTICLE=H-, MASS=0.9, CHARGE=-1, PC=2;",
            None,
            2.1931712199461311,
            0.9,
            -1.0,
        ),
        (
            "BEAM, PARTICLE=PB208+82, MASS=193, CHARGE=82, GAMMA=3;",
            None,
            579.0,
            193.0,
            82.0,
        ),
        ("BEAM, PARTICLE=H- ENERGY=2, MASS=0.9, CHARGE=-1;", None, 2.0, 0.9, -1.0),
        ("BEAM, MASS=2, BRHO=10;", None, 3.60382460552233, 2.0, 1.0),
        # the default particle is a positron at 1 GeV
        ("BEAM;", None, 1.0, 0.00051099895, 1.0),
        ("BEAM, PC=300;", None, 300.0000000004352, 0.00051099895, 1.0),
        ("BEAM, GAMMA=5;", None, 0.00255499475, 0.00051099895, 1.0),
        # a new sequence beam starts from the defaults
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2; BEAM, SEQUENCE=S, PC=1;",
            "s",
            1.0000001305599548,
            0.00051099895,
            1.0,
        ),
        (
            "BEAM, PARTICLE=ELECTRON, ENERGY=0.1; BEAM, SEQUENCE=S, PARTICLE=PROTON;",
            "s",
            1.0,
            0.93827208816,
            1.0,
        ),
        # RESBEAM restores the defaults
        (
            "BEAM, PARTICLE=PROTON, PC=2; RESBEAM; BEAM, PC=1;",
            None,
            1.0000001305599548,
            0.00051099895,
            1.0,
        ),
        (
            "BEAM, PARTICLE=ELECTRON, PC=0.01; RESBEAM; BEAM, PARTICLE=PROTON;",
            None,
            1.0,
            0.93827208816,
            1.0,
        ),
        # a label, commas and flags do not change the beam
        (
            "b1: beam, particle=proton, pc=2;",
            None,
            2.209152441869987,
            0.93827208816,
            1.0,
        ),
        ("beam particle=proton, pc=2;", None, 2.209152441869987, 0.93827208816, 1.0),
        ("BEAM, PARTICLE=PROTON PC=2;", None, 2.209152441869987, 0.93827208816, 1.0),
        (
            "BEAM, PARTICLE=PROTON, -RADIATE, PC=2;",
            None,
            2.209152441869987,
            0.93827208816,
            1.0,
        ),
        # a quoted name is case-sensitive
        (
            "BEAM, PARTICLE='proton', PC=2;",
            None,
            2.209152441869987,
            0.93827208816,
            1.0,
        ),
        ("BEAM, PARTICLE=POSMUON, GAMMA=5;", None, 0.5282918775, 0.1056583755, 1.0),
        ("BEAM, PARTICLE=NEGMUON, GAMMA=5;", None, 0.5282918775, 0.1056583755, -1.0),
    ],
)
def test_beam_matches_madx(script, sequence, energy, mass, charge):
    """Total energy, mass and charge as printed by MAD-X 5.09.03 for this input."""
    parser = MADXParser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MADXInputWarning)
        parser.parse_string(f"{script}\n")
    beam = parser.context.beams[sequence]
    assert beam.energy == pytest.approx(energy, rel=1e-14)
    assert beam.mass == mass
    assert beam.charge == charge


@pytest.mark.parametrize(
    "script",
    [
        "BEAM, PARTICLE=PROTON, ENERGY=0.5;",
        "BEAM, MASS=2, CHARGE=1;",
        "BEAM, PARTICLE=ION, MASS=200, CHARGE=82;",
        "BEAM, PARTICLE=ELECTRON, ENERGY=1; BEAM, MASS=2;",
        # a reference in the command reads the previous mass, of the electron
        "BEAM, PARTICLE=ELECTRON, ENERGY=2; BEAM, PARTICLE=PROTON, ENERGY=4*BEAM->MASS;",
        "BEAM, PARTICLE=PROTON, ENERGY=2; BEAM->MASS=2; BEAM, EX=1;",
    ],
)
def test_beam_energy_must_exceed_mass(script):
    """MAD-X stops on a stated, kept or default energy not above the mass."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MADXInputWarning)
        with pytest.raises(MADXInputError, match="must exceed the mass"):
            MADXParser().parse_string(f"{script}\n")


def test_beam_number_without_a_comma_raises():
    """MAD-X stops on this, which would otherwise drop what follows the number."""
    with pytest.raises(MADXInputError, match="comma"):
        MADXParser().parse_string("BEAM, PARTICLE=PROTON, ENERGY=2 PC=3;\n")


@pytest.mark.parametrize(
    "script",
    [
        "BEAM, PARTICLE=PROTN, PC=2;",
        "BEAM, PARTICLE=FOO, GAMMA=5;",
        "BEAM, PARTICLE=FOO, CHARGE=1, ENERGY=2;",
        'BEAM, PARTICLE="PROTON", PC=2;',
        "BEAM, PARTICLE=FOO, MASS=2, CHARGE=3, PC=1; BEAM, PARTICLE=FOO, PC=1;",
    ],
)
def test_beam_unknown_particle_without_a_mass_raises(script):
    """Stricter than MAD-X, where a typo would turn into a particle of electron mass."""
    with pytest.raises(MADXInputError, match="unknown and has no MASS"):
        MADXParser().parse_string(f"{script}\n")


@pytest.mark.parametrize(
    "script, name, is_nan",
    [
        ("BEAM, PARTICLE=PROTON, PC=sqrt(-1);", "energy", True),
        ("BEAM, PARTICLE=PROTON, GAMMA=asin(2);", "energy", True),
        ("BEAM, PARTICLE=PROTON, BETA=-2;", "energy", True),
        ("BEAM, MASS=sqrt(-1), CHARGE=1, ENERGY=3;", "gamma", True),
        ("BEAM, MASS=2, CHARGE=sqrt(-1), ENERGY=3;", "brho", True),
        ("BEAM, PARTICLE=PROTON, BETA=-1;", "energy", False),
        ("BEAM, PARTICLE=PROTON, BETA=0, EXN=1e-6;", "ex", False),
        ("BEAM, PARTICLE=PROTON, BETA=0, BCURRENT=1e-3, FREQ0=1;", "npart", False),
    ],
)
def test_beam_non_finite_values_match_madx(script, name, is_nan):
    """MAD-X 5.09.03 continues with NaN or infinity, which the parser also keeps."""
    parser = MADXParser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MADXInputWarning)
        parser.parse_string(f"{script}\n")
    value = parser.context.get_beam_attribute("beam", name)
    assert math.isnan(value) if is_nan else value == math.inf


@pytest.mark.parametrize(
    "script, values",
    [
        (
            "BEAM, MASS=-1, CHARGE=1;",
            {"energy": 1.0, "mass": -1.0, "gamma": -1.0, "pc": 0.0},
        ),
        (
            "BEAM, MASS=-1, CHARGE=1, ENERGY=0.5;",
            {"energy": 0.5, "mass": -1.0, "gamma": -0.5},
        ),
        (
            "BEAM, PARTICLE=ION, MASS=-1, CHARGE=1, ENERGY=2;",
            {"energy": 2.0, "mass": -1.0, "gamma": -2.0},
        ),
    ],
)
def test_beam_negative_mass_matches_madx(script, values):
    """MAD-X 5.09.03 reads a negative mass."""
    parser = MADXParser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MADXInputWarning)
        parser.parse_string(f"{script}\n")
    assert not parser.context.beam.is_physical
    for name, value in values.items():
        assert parser.context.get_beam_attribute("beam", name) == value, name


@pytest.mark.parametrize(
    "expression, expected",
    [
        ("sqrt(-1)", math.nan),
        ("log(-1)", math.nan),
        ("log10(-1)", math.nan),
        ("log(0)", -math.inf),
        ("log10(0)", -math.inf),
        ("asin(2)", math.nan),
        ("acos(2)", math.nan),
        ("exp(1000)", math.inf),
        ("sinh(1000)", math.inf),
        ("sinh(-1000)", -math.inf),
        ("cosh(1000)", math.inf),
        ("cosh(-1000)", math.inf),
        ("floor(1e400)", math.inf),
        ("round(1e400)", math.inf),
        ("0^(-1)", math.inf),
        ("(-2)^0.5", math.nan),
        ("(-8)^(1/3)", math.nan),
        ("10^400", math.inf),
        ("(-10)^401", -math.inf),
        ("(-10)^400", math.inf),
        ("0^0", 1.0),
        ("frac(2.5)", 0.5),
        ("frac(-2.5)", -0.5),
        ("frac(1e400)", 0.0),
        ("tgauss(sqrt(-1))", 0.0),
        ("tgauss(0)", 0.0),
        ("tgauss(-1)", 0.0),
    ],
)
def test_math_edges_match_madx(expression, expected):
    """MAD-X 5.09.03 continues with what its C math functions return."""
    parser = MADXParser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MADXInputWarning)
        parser.parse_string(f"X = {expression};\n")
    assert parser.getVariable("x") == pytest.approx(expected, nan_ok=True)


def test_beam_known_particle_ignores_invalid_mass():
    """MAD-X does not read the MASS and CHARGE of a known particle."""
    parser = MADXParser()
    parser.parse_string("BEAM, PARTICLE=PROTON, MASS=sqrt(-1), PC=2;\n")
    assert parser.context.beam.mass == 0.93827208816


def test_beam_messages_name_the_line_of_the_command():
    """Errors and warnings point at the BEAM command, not the next statement."""
    with pytest.raises(MADXInputError) as error:
        MADXParser().parse_string(
            "A = 1;\n\nBEAM, PARTICLE=PROTON, PC=0;\n\n\nQ: QUADRUPOLE, L=1;\n"
        )
    assert error.value.line_number == 3

    with pytest.warns(MADXInputWarning, match=r":2: BEAM sets the energy"):
        MADXParser().parse_string(
            "\nBEAM, PARTICLE=PROTON, ENERGY=3, PC=2;\n\n\nQ: QUADRUPOLE, L=1;\n"
        )


def test_beam_defined():
    """MAD-X stops on a sequence without BEAM, so the parser tracks if there was one."""
    parser = MADXParser()
    assert not parser.context.beam.defined
    parser.parse_string("BEAM;\n")
    assert parser.context.beam.defined


def test_resbeam():
    """RESBEAM restores all parameters of the beam it names, and only those."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=2, BV=-1, FREQ0=5;
BEAM, SEQUENCE=Z, PARTICLE=PROTON, ENERGY=3;
RESBEAM, SEQUENCE=Z;
"""
    )
    z = parser.context.beams["z"]
    assert (z.particle, z.mass, z.energy) == ("positron", 0.51099895000e-3, 1.0)
    assert not z.energy_is_explicit
    default = parser.context.beams[None]
    assert (default.energy, default.bv, default.freq0) == (2.0, -1.0, 5.0e6)

    parser.parse_string("RESBEAM;\n")
    default = parser.context.beams[None]
    assert (default.energy, default.bv, default.freq0) == (1.0, 1.0, 0.0)


def test_beam_attribute_references():
    """`beam->attr` and `beam%sequence->attr` give the values MAD-X 5.09.03 prints."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, PC=2, FREQ0=5;
BEAM, SEQUENCE=Z, PARTICLE=ELECTRON, ENERGY=3;
MOMENTUM = BEAM->PC;
TOTAL = BEAM->ENERGY;
RIGIDITY = BEAM->BRHO;
FREQUENCY = BEAM->FREQ0;
TOTAL_Z = BEAM%Z->ENERGY;
TOTAL_NONE = BEAM%NOSUCH->ENERGY;
"""
    )
    get = parser.getVariable
    assert get("momentum") == pytest.approx(2.0, rel=1e-14)
    assert get("total") == pytest.approx(2.2091524418699868, rel=1e-14)
    assert get("rigidity") == pytest.approx(6.67128190396304, rel=1e-14)
    assert get("frequency") == pytest.approx(5.0, rel=1e-14)
    assert get("total_z") == 3.0
    assert get("total_none") == pytest.approx(2.2091524418699868, rel=1e-14)


@pytest.mark.parametrize(
    "script, sequence, values",
    [
        # all defaults: no command, and an empty one
        (
            "",
            None,
            {
                "pc": 0.999999869440028,
                "gamma": 1956.9511835591834,
                "beta": 0.999999869440028,
                "brho": 3.335640516480331,
                "ex": 1.0,
                "exn": 1956.9509280596917,
                "ey": 1.0,
                "eyn": 1956.9509280596917,
                "et": 0.001,
                "sigt": 1.0,
                "sige": 0.001,
                "kbunch": 1.0,
                "npart": 1.0,
                "bcurrent": 0.0,
                "freq0": 0.0,
                "circ": 0.0,
                "dtbyds": 0.0,
                "deltap": 0.0,
                "alfa": 2.611199269011025e-07,
                "u0": 0.0,
                "qs": 0.0,
                "arad": 2.8179403246496351e-15,
                "bv": 1.0,
                "bunched": 1.0,
                "radiate": 0.0,
                "n1min": -1.0,
                "pdamp": 0.0,
            },
        ),
        # stated values, of which a BEAM command ignores DELTAP to ARAD
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, EX=2e-6, EY=3e-6, ET=4e-3, SIGT=0.5, "
            "SIGE=2e-4, KBUNCH=7, NPART=1e10, DELTAP=2, U0=3, QS=4, N1MIN=5, ALFA=6, "
            "ARAD=7;",
            None,
            {
                "ex": 2e-06,
                "exn": 3.764903254435854e-06,
                "ey": 3e-06,
                "eyn": 5.647354881653782e-06,
                "et": 0.004,
                "sigt": 0.5,
                "sige": 0.0002,
                "kbunch": 7.0,
                "npart": 1e10,
                "bcurrent": 0.0,
                "deltap": 0.0,
                "u0": 0.0,
                "qs": 0.0,
                "n1min": -1.0,
                "alfa": 0.22008862785503169,
                "arad": 1.5346982663445394e-18,
            },
        ),
        # normalized emittances, of which a later command keeps the unnormalized
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, EXN=1e-6, EYN=2e-6;",
            None,
            {
                "ex": 5.312221496378628e-07,
                "exn": 1e-06,
                "ey": 1.0624442992757256e-06,
                "eyn": 2e-06,
            },
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, EXN=1e-6; BEAM, ENERGY=3;",
            None,
            {"ex": 5.312221496378628e-07, "exn": 1.6133028901394618e-06},
        ),
        # CIRC sets FREQ0 and takes precedence over it; a later energy keeps both
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, CIRC=100;",
            None,
            {"freq0": 2.647544871519702, "circ": 100.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, CIRC=100, FREQ0=5;",
            None,
            {"freq0": 2.647544871519702, "circ": 100.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, FREQ0=5;",
            None,
            {"freq0": 5.0, "circ": 52.95089743039404},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, CIRC=100; BEAM, ENERGY=3;",
            None,
            {"freq0": 2.647544871519702, "circ": 100.0},
        ),
        # BCURRENT and NPART follow from each other; each command clears BCURRENT
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, BCURRENT=1e-3, FREQ0=1;",
            None,
            {"npart": 7067518919.850619, "bcurrent": 0.001, "circ": 264.7544871519702},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, BCURRENT=1e-3, FREQ0=1; BEAM, ENERGY=3;",
            None,
            {"npart": 7067518919.850619, "bcurrent": 0.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2, NPART=5, FREQ0=2;",
            None,
            {
                "npart": 5.0,
                "bcurrent": 1.4149236971850033e-12,
                "circ": 132.3772435759851,
            },
        ),
        # a stated energy parameter is kept as stated
        (
            "BEAM, PARTICLE=PROTON, PC=1e-9;",
            None,
            {"pc": 1e-09, "beta": 1.0657889247894517e-09, "gamma": 1.0},
        ),
        ("BEAM, PARTICLE=PROTON, PC=1e-6;", None, {"pc": 1e-06}),
        (
            "BEAM, PARTICLE=PROTON, BETA=-0.5;",
            None,
            {
                "beta": -0.5,
                "pc": 0.5417116426722884,
                "energy": 1.0834232853445767,
            },
        ),
        ("BEAM, PARTICLE=PROTON, BETA=1e-9;", None, {"beta": 1e-09}),
        (
            "BEAM, PARTICLE=PROTON, BRHO=2;",
            None,
            {"brho": 2.0, "pc": 0.5995849160000001, "energy": 1.1134884745316647},
        ),
        # a parameter without a value takes its default
        (
            "BEAM, PARTICLE=PROTON, ENERGY=5; BEAM, BETA;",
            None,
            {"energy": 0.93827208816, "beta": 0.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=5; BEAM, -BRHO;",
            None,
            {"energy": 0.93827208816, "brho": 0.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=5; BEAM, PARTICLE;",
            None,
            {"mass": 0.00051099895, "charge": 1.0, "energy": 5.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=5; BEAM, MASS;",
            None,
            {"mass": 0.93827208816, "charge": 1.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=5, EX=2; BEAM, EX;",
            None,
            {"ex": 0.0, "exn": 0.0},
        ),
        (
            "BEAM, SEQUENCE, PARTICLE=PROTON, ENERGY=3;",
            None,
            {"energy": 3.0, "mass": 0.93827208816},
        ),
        ("BEAM, RADIATE;", None, {"radiate": 1.0, "bunched": 1.0}),
        ("BEAM, -RADIATE, -BUNCHED;", None, {"radiate": 0.0, "bunched": 0.0}),
        ("BEAM, RADIATE=TRUE, BUNCHED=FALSE;", None, {"radiate": 1.0, "bunched": 0.0}),
        # a quoted sequence name keeps its case, so `a` keeps the default beam
        (
            'BEAM, PARTICLE=ELECTRON, ENERGY=2; BEAM, SEQUENCE="A", PARTICLE=PROTON, ENERGY=5;',
            "a",
            {"energy": 2.0, "mass": 0.00051099895},
        ),
        (
            'BEAM, SEQUENCE=X, PARTICLE=ELECTRON, ENERGY=3; RESBEAM, SEQUENCE="X";',
            "x",
            {"energy": 3.0, "mass": 0.00051099895},
        ),
        # a reference in a BEAM command reads the beam before this command
        (
            "BEAM, PARTICLE=PROTON, ENERGY=3; BEAM, MASS=2, CHARGE=1, ENERGY=4*BEAM->MASS;",
            None,
            {"energy": 3.7530883526399998, "mass": 2.0, "charge": 1.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=3; BEAM, MASS=2, CHARGE=1, ENERGY:=4*BEAM->MASS;",
            None,
            {"energy": 3.7530883526399998, "mass": 2.0, "charge": 1.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=3; BEAM, PC=BEAM->ENERGY;",
            None,
            {"energy": 3.1433031211482176, "pc": 3.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=3, EX=2; "
            "BEAM, ENERGY=5, EY=BEAM->ENERGY, ET=BEAM->EXN;",
            None,
            {
                "energy": 5.0,
                "ey": 3.0,
                "et": 6.073929301476829,
                "exn": 10.46855306240483,
            },
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=3; "
            "BEAM, ENERGY=5, CIRC=BEAM->ENERGY, NPART=BEAM->GAMMA;",
            None,
            {"energy": 5.0, "circ": 3.0, "npart": 3.1973667743683549},
        ),
        (
            "BEAM, PARTICLE=ION, MASS=2*BEAM->MASS, CHARGE=3, ENERGY=9;",
            None,
            {"mass": 0.0010219979, "charge": 3.0, "energy": 9.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=3; "
            "BEAM, MASS=BEAM->ENERGY, CHARGE=BEAM->MASS, ENERGY=9;",
            None,
            {"mass": 3.0, "charge": 0.93827208816, "energy": 9.0},
        ),
        (
            "BEAM, SEQUENCE=S, PARTICLE=PROTON, ENERGY=3; "
            "BEAM, SEQUENCE=S, PC=BEAM%S->ENERGY;",
            "s",
            {"energy": 3.1433031211482176, "pc": 3.0},
        ),
        # ... which falls back to the default beam for a new sequence
        (
            "BEAM, PARTICLE=PROTON, ENERGY=3; "
            "BEAM, SEQUENCE=S, PARTICLE=PROTON, ENERGY=2*BEAM%S->ENERGY;",
            "s",
            {"energy": 6.0},
        ),
        # an assignment stores only this value, until the next BEAM command
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2; BEAM->ENERGY=7;",
            None,
            {
                "energy": 7.0,
                "pc": 1.7662518191299543,
                "gamma": 2.1315778495789033,
            },
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2; BEAM->ENERGY=7; BEAM, EX=1;",
            None,
            {"energy": 7.0, "pc": 6.9368325256257899},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2; X=3; BEAM->EX:=X; X=5;",
            None,
            {"ex": 3.0},
        ),
        (
            "BEAM, SEQUENCE=S, PARTICLE=PROTON, ENERGY=2; BEAM%S->ENERGY=7;",
            "s",
            {"energy": 7.0},
        ),
        (
            "BEAM, SEQUENCE=S, PARTICLE=PROTON, ENERGY=2; BEAM%S->ENERGY=7;",
            None,
            {"energy": 1.0},
        ),
        (
            "BEAM, PARTICLE=PROTON, ENERGY=2; BEAM%NOSUCH->ENERGY=7;",
            None,
            {"energy": 7.0},
        ),
        ("BEAM->ENERGY=7;", None, {"energy": 7.0}),
        ("BEAM, PARTICLE=PROTON, ENERGY=2; BEAM->FREQ0=7;", None, {"freq0": 7.0}),
        ("BEAM, PARTICLE=PROTON, ENERGY=2; BEAM->FOO=7;", None, {"foo": 0.0}),
        # a label does not change a command, but BEAM and RESBEAM are protected
        (
            "BEAM, PARTICLE=PROTON, PC=2; r: resbeam;",
            None,
            {"energy": 1.0, "mass": 0.00051099895},
        ),
        (
            "BEAM, PARTICLE=PROTON, PC=2; beam: resbeam;",
            None,
            {"energy": 2.209152441869987, "mass": 0.93827208816},
        ),
        (
            "beam: beam, particle=proton, pc=2;",
            None,
            {"energy": 1.0, "mass": 0.00051099895},
        ),
        (
            "resbeam: beam, particle=proton, pc=2;",
            None,
            {"energy": 1.0, "mass": 0.00051099895},
        ),
    ],
)
def test_beam_references_match_madx(script, sequence, values):
    """`beam->name` gives what MAD-X 5.09.03 prints with `value` for this input."""
    parser = MADXParser()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MADXInputWarning)
        parser.parse_string(f"{script}\n")
    prefix = f"beam%{sequence}" if sequence else "beam"
    for name, value in values.items():
        got = parser.context.get_beam_attribute(prefix, name)
        assert got == pytest.approx(value, rel=1e-13, abs=0.0), name


@pytest.mark.parametrize(
    "script, match",
    [
        ("BEAM, ENERGIE=3;", "Illegal keyword: energie"),
        ("BEAM, E=2;", "Illegal keyword: e"),
        ("BEAM, PART=PROTON, ENERGY=2;", "Illegal keyword: part"),
        ("BEAM, PARTICLE=PROTON, ENERGY=2, FOO=3;", "Illegal keyword: foo"),
        ("RESBEAM, ENERGY=3;", "Illegal keyword: energy"),
    ],
)
def test_beam_illegal_keyword_raises(script, match):
    """MAD-X stops on an unknown parameter, which must not fall back to 1 GeV."""
    with pytest.raises(MADXInputError, match=match):
        MADXParser().parse_string(f"{script}\n")


def test_beam_illegal_keyword_after_a_string_raises():
    """Stricter than MAD-X, which reads the unknown particle `proton,energie=3`."""
    with pytest.raises(MADXInputError, match="Illegal keyword: energie"):
        MADXParser().parse_string("BEAM, PARTICLE=PROTON, ENERGIE=3;\n")


@pytest.mark.parametrize(
    "param, match",
    [
        ("ENERGY", "must exceed the mass"),
        ("-ENERGY", "must exceed the mass"),
        ("PC", "PC must be > 0"),
        ("-PC", "PC must be > 0"),
        ("GAMMA", "GAMMA must be > 1"),
    ],
)
def test_beam_parameter_without_a_value_stops(param, match):
    """MAD-X reads these as zero, and stops as for a stated zero."""
    with pytest.raises(MADXInputError, match=match):
        MADXParser().parse_string(f"BEAM, PARTICLE=PROTON, ENERGY=5; BEAM, {param};\n")


def test_predefined_constants_match_madx():
    """Predefined constants as printed by MAD-X 5.09.03."""
    parser = MADXParser()
    for name, value in {
        "pi": 3.1415926535897931,
        "twopi": 6.2831853071795862,
        "degrad": 57.295779513082323,
        "raddeg": 0.017453292519943295,
        "e": 2.7182818284590451,
        "amu0": 1.2566370614359173e-06,
        "emass": 0.00051099895000000003,
        "mumass": 0.1056583755,
        "nmass": 0.93956542052000003,
        "umass": 0.93149410241999997,
        "pmass": 0.93827208815999996,
        "clight": 299792458,
        "qelect": 1.6021766339999999e-19,
        "hbar": 6.5821195689999996e-25,
        "erad": 2.8179403262e-15,
        "prad": 1.5346982671888944e-18,
    }.items():
        assert parser.getVariable(name) == pytest.approx(value, rel=1e-15), name


@pytest.mark.parametrize(
    "script, name, value",
    [
        ("PI=3;", "pi", math.pi),
        ("PI:=3;", "pi", math.pi),
        ("AMU0=1;", "amu0", 4.0e-7 * math.pi),
        ("CONST EMASS=3;", "emass", 0.51099895000e-3),
        ("CONST C1=2; C1=3;", "c1", 2.0),
    ],
)
def test_redefining_a_constant_is_ignored(script, name, value):
    """MAD-X warns and continues with the constant."""
    parser = MADXParser()
    with pytest.warns(MADXInputWarning, match="redefine constant"):
        parser.parse_string(f"{script}\n")
    assert parser.getVariable(name) == value


def test_beam_assignment_does_not_state_the_energy():
    """What follows from an assigned value is out of date, as in MAD-X."""
    parser = MADXParser()
    parser.parse_string("BEAM, PARTICLE=PROTON, PC=2;\n")
    assert parser.context.beam.energy_is_explicit
    parser.parse_string("BEAM->ENERGY=7;\n")
    assert not parser.context.beam.energy_is_explicit


@pytest.mark.parametrize(
    "statement, skipped",
    [("USE, SEQUENCE=NOPE;", "nope"), ('USE, SEQUENCE="A";', "A")],
)
def test_use_skips_an_unknown_sequence(statement, skipped):
    """MAD-X warns and keeps the current sequence, and a quoted name keeps its case."""
    parser = MADXParser()
    parser.parse_string("D1: DRIFT, L=1;\nA: LINE=(D1);\nBEAM;\nUSE, SEQUENCE=A;\n")
    with pytest.warns(MADXInputWarning, match=f"unknown sequence skipped: {skipped}"):
        parser.parse_string(f"{statement}\n")
    assert parser.sequence["name"] == "a"


def test_lenient_commas():
    """Test that missing commas in element and LINE definitions are tolerated."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=7000.0;
! Missing comma between K1 and K1S attributes
QSS01: QUADRUPOLE, L=0.5, K1=1.0 K1S:=0.0;
D1: DRIFT, L=1.0;
! Missing comma between LINE elements
MYLINE: LINE=(QSS01 D1 QSS01);
USE, SEQUENCE=MYLINE;
"""
    )
    beamline = parser.getBeamline()
    assert len(beamline) == 3
    assert beamline[0]["type"] == "quadrupole"
    assert beamline[1]["type"] == "drift"
    assert beamline[2]["type"] == "quadrupole"


def test_exist_function():
    """Test EXIST() function in IF conditions."""
    parser = MADXParser()
    parser.parse_string(
        """\
myvar = 42;
found = 0;
missing = 0;
IF (EXIST(myvar)) { found = 1; };
IF (EXIST(novar)) { missing = 1; };
"""
    )
    assert parser.context.get_variable("found") == 1.0
    assert parser.context.get_variable("missing") == 0.0


def test_deferred_expressions():
    """Test that := creates deferred expressions re-evaluated on access."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=7000.0;
k = 1.0;
QF: QUADRUPOLE, L=0.5, K1:=k;
QD: QUADRUPOLE, L=0.5, K1:=-k;
MYLINE: LINE=(QF, QD);
! Change k after element definitions
k = 3.0;
USE, SEQUENCE=MYLINE;
"""
    )
    beamline = parser.getBeamline()
    # Deferred := means K1 picks up the updated k=3.0
    assert beamline[0]["k1"] == 3.0
    assert beamline[1]["k1"] == -3.0


def test_predefined_constants():
    """Test that predefined MAD-X constants (pi, clight, etc.) are available.

    Regression test for a bug where constants were dead code after a return statement.
    """
    parser = MADXParser()
    parser.parse_string(
        """\
x = pi;
y = twopi;
z = clight;
m = pmass;
"""
    )
    assert parser.context.get_variable("x") == pytest.approx(math.pi)
    assert parser.context.get_variable("y") == pytest.approx(2 * math.pi)
    assert parser.context.get_variable("z") == pytest.approx(2.99792458e8)
    assert parser.context.get_variable("m") == pytest.approx(0.93827208816)


def test_get_beamline_without_use_infers_unique_root_line():
    """Declaration-only files may omit USE; infer the unique top-level line."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=1.0;
D1: DRIFT, L=0.5;
CELL: LINE=(D1, D1);
BOOSTER: LINE=(CELL, CELL);
"""
    )

    with pytest.warns(
        MADXInputWarning, match="inferring unique top-level line 'booster'"
    ):
        beamline = parser.getBeamline()

    assert [elem["type"] for elem in beamline] == ["drift"] * 4
    assert parser.sequence["name"] == "booster"


def test_get_beamline_without_use_requires_explicit_choice_when_roots_ambiguous():
    """Do not guess when multiple top-level root lines exist."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=PROTON, ENERGY=1.0;
D1: DRIFT, L=0.5;
A: LINE=(D1);
B: LINE=(D1, D1);
"""
    )

    with pytest.raises(MADXInputError, match="multiple top-level lines"):
        parser.getBeamline()

    beamline = parser.getBeamline(line_name="b")
    assert [elem["type"] for elem in beamline] == ["drift", "drift"]


# ---------------------------------------------------------------------------
# Line continuation ('&' at end-of-line)
# ---------------------------------------------------------------------------


def test_continuation_in_expression():
    """A long arithmetic expression split across lines with '&' parses
    identically to the single-line form."""
    parser = MADXParser()
    parser.parse_string(
        """\
LEN = 1 + 2 &
    + 3 &
    + 4;
"""
    )
    assert parser.getVariable("len") == pytest.approx(10.0)


def test_continuation_in_line_definition():
    """A LINE = (...) definition split across multiple lines with '&'
    expands to the same sequence as the one-line equivalent."""
    parser = MADXParser()
    parser.parse_string(
        """\
BEAM, PARTICLE=ELECTRON, ENERGY=2.0;
D1: DRIFT, L=0.1;
QF: QUADRUPOLE, L=0.1, K1=1.0;
QD: QUADRUPOLE, L=0.1, K1=-1.0;
FODO: LINE=(D1, QF, &
            D1, QD, &
            D1, QF);
USE, SEQUENCE=FODO;
"""
    )
    beamline = parser.getBeamline()
    assert len(beamline) == 6
    assert [e["type"] for e in beamline] == [
        "drift",
        "quadrupole",
        "drift",
        "quadrupole",
        "drift",
        "quadrupole",
    ]
    assert beamline[1]["k1"] == pytest.approx(1.0)
    assert beamline[3]["k1"] == pytest.approx(-1.0)


def test_continuation_with_trailing_whitespace():
    """Trailing spaces/tabs between the value and '&' are allowed."""
    parser = MADXParser()
    parser.parse_string("X = 1 &   \t  \n  + 2;\n")
    assert parser.getVariable("x") == pytest.approx(3.0)


def test_continuation_with_trailing_comment():
    """A trailing line-comment ('!' or '//') after '&' is allowed and
    ignored, matching MAD-X's read-then-strip behavior."""
    parser = MADXParser()
    parser.parse_string(
        """\
Y = 10 & ! this is a continuation
    + 5 & // and another one
    + 1;
"""
    )
    assert parser.getVariable("y") == pytest.approx(16.0)


def test_continuation_chained():
    """Five or more chained '&' continuations join into one statement."""
    parser = MADXParser()
    parser.parse_string(
        """\
Z = 1 &
  + 1 &
  + 1 &
  + 1 &
  + 1 &
  + 1;
"""
    )
    assert parser.getVariable("z") == pytest.approx(6.0)


def test_continuation_midline_raises():
    """A '&' not at end-of-line (and not part of '&&') must raise
    MADXInputError with a helpful message."""
    parser = MADXParser()
    with pytest.raises(MADXInputError, match="line-continuation"):
        parser.parse_string("A = 1 & 2;\n")


def tight_approx(expected):
    """Agree to rounding: ``rel`` alone keeps the default ``abs`` of 1e-12."""
    return pytest.approx(expected, rel=1e-15, abs=0.0)


@pytest.mark.parametrize(
    "expression, expected",
    [
        ("(2^3)^2", 64.0),
        ("2^(3^2)", 512.0),
        ("2^(-1)", 0.5),
        ("-2^2", -4.0),
        ("-2^2*3", -12.0),
        ("-2^2+1", -3.0),
        ("(-2^2)", -4.0),
        ("2*(-3^2)", -18.0),
        ("sqrt(-2^2+8)", 2.0),
        ("-2^0.5", -1.4142135623730951),
        ("(-2)^2", 4.0),
        ("-(2^2)", -4.0),
        ("1-2^2", -3.0),
        ("-2*3", -6.0),
        ("3*(-2)", -6.0),
        ("-(-2)", 2.0),
        ("2*exp(-1)", 0.73575888234288467),
    ],
)
def test_power_as_in_madx(expression, expected):
    """``^`` binds tighter than a sign, as in MAD-X.

    A sign at the start of an expression and after ``(`` applies to the whole
    power. See polish_expr and polish_value in src/mad_eval.c and pre_split in
    src/mad_parse.c. The expected values are from MAD-X 5.09.03.
    """
    parser = MADXParser()
    parser.parse_string(f"X = {expression};\n")
    assert parser.getVariable("x") == tight_approx(expected)


@pytest.mark.parametrize(
    "expression, pattern",
    [
        ("3*-2", "*-"),
        ("3/-2", "/-"),
        ("3*+2", "*+"),
        ("2--3", "--"),
        ("2+-3", "+-"),
        ("3-+2", "-+"),
        ("2^-1", "^-"),
        ("-2^-2", "^-"),
        ("3 * - 2", "*-"),
        ("10 + 3*-2", "*-"),
        ("(3*-2)", "*-"),
        ("(3--2)", "--"),
        ("(2^-1)", "^-"),
        ("sqrt(2^-2)", "^-"),
    ],
)
def test_sign_after_operator_raises(expression, pattern):
    """A sign directly after an operator raises, as MAD-X does not evaluate it as written.

    MAD-X rejects these outside of parentheses (loc_expr, src/mad_expr.c).
    Inside of them, it applies an operator that lacks its left operand to 0
    (polish_value, src/mad_eval.c): in MAD-X 5.09.03, (3*-2) is -2, (3--2) is
    -5, (2^-1) is -1 and sqrt(2^-2) is NaN.
    """
    parser = MADXParser()
    with pytest.raises(
        MADXInputError,
        match=re.escape(f'a sign directly after an operator, "{pattern}"'),
    ):
        parser.parse_string(f"X = {expression};\n")


@pytest.mark.parametrize(
    "expression, pattern",
    [
        ("--2", "--"),
        ("+-2", "+-"),
        ("(--2)", "--"),
        ("1 + (-+2)", "-+"),
    ],
)
def test_two_leading_signs_raise(expression, pattern):
    """Two signs in a row raise. In MAD-X 5.09.03, (--2) is -2."""
    parser = MADXParser()
    with pytest.raises(
        MADXInputError, match=re.escape(f'two signs in a row, "{pattern}"')
    ):
        parser.parse_string(f"X = {expression};\n")


@pytest.mark.parametrize(
    "expression",
    [
        "2^3^2",
        "2^2^3^2",
        "3*2^3^2",
        "2^3^2*3",
        "-2^3^2",
        "(2^3)^2^2",
        "sqrt(4)^2^2",
    ],
)
def test_chained_power_raises(expression):
    """A chained power raises: MAD-X evaluates it from the left, 2^3^2 = 64.

    See polish_expr in src/mad_eval.c; the usual reading is 2^(3^2) = 512.
    """
    parser = MADXParser()
    with pytest.raises(MADXInputError, match="a chained power"):
        parser.parse_string(f"X = {expression};\n")


@pytest.mark.parametrize(
    "statements",
    [
        "A = 3; B = 2; X := A*-B;",
        "CONST X = 2*-1;",
        "INT N = 2*-1;",
        "A = 3; B = 2; Q: QUADRUPOLE, L=1, K1=(A*-B);",
        "A = 3; B = 2; Q: QUADRUPOLE, L=1, K1:=A^B^2;",
        "M: MULTIPOLE, KNL={0, 2*-1};",
        "S: SEQUENCE, L=10; Q: QUADRUPOLE, L=1, AT=5*-1; ENDSEQUENCE;",
        "BEAM, PARTICLE=PROTON, ENERGY=2*-1;",
        "BEAM, PARTICLE=PROTON, ENERGY=2; BEAM->ENERGY = 2*-1;",
        "X = 1; IF (X > 2*-1) { Y = 1; }",
        "X = 1; WHILE (X < 2^2^2) { X = X + 1; }",
        "X = MAX(1, 2*-1);",
        'X := 1; ADD2EXPR, VAR=X, EXPR="2*-1";',
    ],
)
def test_unsupported_expression_contexts(statements):
    """Every kind of statement raises for these forms, also deferred ones."""
    parser = MADXParser()
    with pytest.raises(MADXInputError, match="is not supported"):
        parser.parse_string(statements + "\n")


@pytest.mark.parametrize(
    "expression, expected",
    [
        ("-2", -2.0),
        ("+2", 2.0),
        ("3*(+2)", 6.0),
        ("max(-1, -2)", -1.0),  # ImpactX extension: MAD-X has no max
        ("1e-3*2", 0.002),
        ("1.d-3*2", 0.002),
        ("1E+3*2", 2000.0),
    ],
)
def test_leading_sign_accepted(expression, expected):
    """A sign at the start, after ``(`` or after ``,`` stays accepted, as do exponents."""
    parser = MADXParser()
    parser.parse_string(f"X = {expression};\n")
    assert parser.getVariable("x") == tight_approx(expected)


def test_unsupported_expression_location(tmp_madx):
    """The error names the file and line, also in a file read with CALL, in a
    macro defined there and in a quoted ADD2EXPR expression."""
    main = tmp_madx("main.madx", "A = 1;\nB = 2^-1;\n")
    with pytest.raises(MADXInputError) as excinfo:
        MADXParser().parse(main)
    assert excinfo.value.file == main
    assert excinfo.value.line_number == 2
    assert str(excinfo.value).startswith(f"{main}:2: ")

    main = tmp_madx("call.madx", "A = 1;\nCALL, FILE='strengths.str';\nB = 2;\n")
    strengths = tmp_madx("strengths.str", "KQ = 0.5;\n\nKQD := KQ^2^2;\n")
    with pytest.raises(MADXInputError) as excinfo:
        MADXParser().parse(main)
    assert excinfo.value.file == strengths
    assert excinfo.value.line_number == 3

    for expr in ("2*-1", "2^2^2"):
        main = tmp_madx(
            "add2expr.madx", f'X := 1;\n\nADD2EXPR, VAR=X, EXPR="{expr}";\n'
        )
        with pytest.raises(MADXInputError, match="is not supported") as excinfo:
            MADXParser().parse(main)
        assert excinfo.value.file == main
        assert excinfo.value.line_number == 3

    # the macro is defined in macros.str and executed from exec.madx
    main = tmp_madx(
        "exec.madx", "CALL, FILE='macros.str';\nB = 2;\nKQ = 3;\nEXEC, SETK(KQ);\n"
    )
    macros = tmp_madx("macros.str", "A = 1;\nSETK(V): MACRO = { K = V*-1; };\n")
    with pytest.raises(MADXInputError, match="is not supported") as excinfo:
        MADXParser().parse(main)
    assert excinfo.value.file == macros
    assert excinfo.value.line_number == 2


def test_incomplete_expression_location(tmp_madx):
    """An expression that ends early names its line, not line 0."""
    main = tmp_madx("main.madx", "A = 1;\nB = 2+;\n")
    with pytest.raises(MADXInputError, match="Unexpected token") as excinfo:
        MADXParser().parse(main)
    assert excinfo.value.file == main
    assert excinfo.value.line_number == 2


def test_elseif_after_taken_branch_not_evaluated():
    """As in MAD-X, no ELSEIF condition is evaluated after a taken branch.

    MAD-X 5.09.03 gives Y = 1 here, without a warning.
    """
    parser = MADXParser()
    parser.parse_string(
        "X = 1;\nIF (X > 0) { Y = 1; } ELSEIF (X > 2*-1) { Y = 2; } ELSE { Y = 3; }\n"
    )
    assert parser.getVariable("y") == 1.0

    # a reached ELSEIF condition is evaluated
    with pytest.raises(MADXInputError, match="is not supported"):
        MADXParser().parse_string(
            "X = -5;\nIF (X > 0) { Y = 1; } ELSEIF (X > 2*-1) { Y = 2; }\n"
        )


def test_power_parentheses_deferred():
    """Nested powers through variables, a deferred expression and an element."""
    parser = MADXParser()
    parser.parse_string(
        """\
A = 2;
B = 3;
NOW = (A^B)^2;
LATER := (A^B)^2;
QF: QUADRUPOLE, L=1.0, K1:=(A^B)^2;
FODO: LINE=(QF);
USE, SEQUENCE=FODO;
A = 3;
"""
    )
    assert parser.getVariable("now") == tight_approx(64.0)
    # the deferred expressions see the updated A: (3^3)^2
    assert parser.getVariable("later") == tight_approx(729.0)
    assert parser.getBeamline()[0]["k1"] == tight_approx(729.0)


def test_power_sign_deferred():
    """A deferred ``-a^2`` follows ``a`` and stays negative, as in MAD-X."""
    parser = MADXParser()
    parser.parse_string(
        """\
A = 2;
X := -A^2;
AT2 = X;
A = 3;
AT3 = X;
"""
    )
    assert parser.getVariable("at2") == tight_approx(-4.0)
    assert parser.getVariable("at3") == tight_approx(-9.0)
    assert parser.getVariable("x") == tight_approx(-9.0)


def test_power_sign_element_parameter():
    """A sign in an element parameter applies to the whole power, as in MAD-X."""
    parser = MADXParser()
    parser.parse_string(
        """\
QF: QUADRUPOLE, L=1.0, K1=-2^2;
FODO: LINE=(QF);
USE, SEQUENCE=FODO;
"""
    )
    assert parser.getBeamline()[0]["k1"] == tight_approx(-4.0)


def test_power_sign_condition():
    """IF and WHILE conditions use the same grammar, as in MAD-X: -2^2 is -4."""
    parser = MADXParser()
    parser.parse_string(
        """\
IF (-2^2 < 0) { Y1 = 1; } ELSE { Y1 = 2; };
IF (-2^2 == -4) { Y2 = 1; } ELSE { Y2 = 2; };
IF (-2^2 > 0) { Y3 = 1; } ELSEIF (-2^2 < 0) { Y3 = 3; } ELSE { Y3 = 2; };
N = 0;
WHILE (-(N+1)^2 > -10) { N = N + 1; };
"""
    )
    assert parser.getVariable("y1") == 1.0
    assert parser.getVariable("y2") == 1.0
    assert parser.getVariable("y3") == 3.0
    assert parser.getVariable("n") == 3.0


def test_power_sign_add2expr():
    """A quoted ``ADD2EXPR`` expression uses the same grammar as other input."""
    parser = MADXParser()
    parser.parse_string(
        """\
A = 2;
X := 1;
ADD2EXPR, VAR=X, EXPR="-A^2";
"""
    )
    assert parser.getVariable("x") == tight_approx(-3.0)
