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

import pytest

from impactx.MADXParser import MADXInputError, MADXParser


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
    with pytest.raises(FileNotFoundError):
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
