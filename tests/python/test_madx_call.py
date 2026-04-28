#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""
Tests for the MAD-X CALL command support.

References:
    - https://github.com/BLAST-ImpactX/impactx/issues/512
    - http://madx.web.cern.ch/madx/webguide/manual.html (Chapter 5.2)
"""

import math

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
