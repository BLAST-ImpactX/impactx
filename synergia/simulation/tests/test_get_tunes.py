#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest


# lattice_fixture is a lattice
@pytest.fixture
def lattice_fixture():
    # booster-like lattice
    simplefodo_madx = """
f: multipole, knl={0, 1/8.1};
d: multipole, knl={0, -1/9.4};
fodo: sequence, l=20, refer=entry;
  f, at=0.0;
  d, at=10.0;
endsequence;

beam, particle=proton, pc=0.75*pmass;
! tunes:  [ 0.22239535 -0.22239535  0.16407172 -0.16407172]
! alphax:  -1.9200614429492773
! betax:  31.10499537577829
! alphay:  0.6737128494450171
! betay:  10.914148161009276
! chromx:  -0.3996685356 from MADX
! chromy:  -0.3533542246 from MADX


"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(simplefodo_madx)
    lattice = reader.get_lattice("fodo")
    lattice.set_all_string_attribute("extractor_type", "libff")

    return lattice


def test_tunes(lattice_fixture):
    tunes = synergia.simulation.Lattice_simulator.calculate_tune_and_cdt(
        lattice_fixture
    )
    print("Tunes: x: ", tunes[0], ", y: ", tunes[1])
    assert tunes[0] == pytest.approx(0.22239535)
    assert tunes[1] == pytest.approx(0.16407172)


def test_chromaticities(lattice_fixture):
    beta = lattice_fixture.get_reference_particle().get_beta()
    chrom = synergia.simulation.Lattice_simulator.get_chromaticities(lattice_fixture)
    print("H chromaticitity: ", chrom.horizontal_chromaticity)
    print("V chromaticitity: ", chrom.vertical_chromaticity)
    # MADX chromaticiticity is with respect to PT
    assert chrom.horizontal_chromaticity == pytest.approx(-0.3996685356 * beta)
    assert chrom.vertical_chromaticity == pytest.approx(-0.3533542246 * beta)


def main():
    lf = lattice_fixture()
    test_tunes(lf)
    test_chromaticities(lf)


if __name__ == "__main__":
    main()
