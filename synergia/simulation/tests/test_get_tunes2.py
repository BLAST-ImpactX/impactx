#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

# This test is similar to test_get_tunes2.py but with the Booster-like lattice.


# lattice_fixture is a lattice
@pytest.fixture
def lattice_fixture():
    # booster-like lattice
    booster_madx = """
ncells=24;
turn_voltage=1.0; ! 1 MV /turn
beam, particle=proton,energy=pmass+0.8;

f: sbend, l=2.0, angle=(pi/(2*ncells)), k1=1/16.2;
d: sbend, l=2.0, angle=(pi/(2*ncells)), k1=-1/16.7;
!f: quadrupole, l=2.0, k1=0.0625;
!d: quadrupole, l=2.0, k1=-0.0625;
rfc: rfcavity, l=0.0, volt=turn_voltage/ncells, harmon=96, lag=(1/120.0);

cell: sequence, l=20.0, refer=centre;
fodo_1: f, at=1.5;
fodo_2: d, at=8.5;
fodo_3: d, at=11.5;
fodo_4: f, at=18.5;
fodo_5: rfc, at=20.0;
endsequence;

booster: sequence, l=480.0, refer=entry;
cell, at=0.0;
cell, at=20.0;
cell, at=40.0;
cell, at=60.0;
cell, at=80.0;
cell, at=100.0;
cell, at=120.0;
cell, at=140.0;
cell, at=160.0;
cell, at=180.0;
cell, at=200.0;
cell, at=220.0;
cell, at=240.0;
cell, at=260.0;
cell, at=280.0;
cell, at=300.0;
cell, at=320.0;
cell, at=340.0;
cell, at=360.0;
cell, at=380.0;
cell, at=400.0;
cell, at=420.0;
cell, at=440.0;
cell, at=460.0;
endsequence;

"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(booster_madx)
    lattice = reader.get_lattice("booster")
    lattice.set_all_string_attribute("extractor_type", "libff")

    return lattice


def test_tunes(lattice_fixture):
    tunes = synergia.simulation.Lattice_simulator.calculate_tune_and_cdt(
        lattice_fixture
    )
    print("Tunes: x: ", tunes[0], ", y: ", tunes[1])
    assert tunes[0] == pytest.approx(0.16533417, rel=1.0e-4)
    assert tunes[1] == pytest.approx(0.436092762, rel=1.0e-4)


def test_chromaticities(lattice_fixture):
    beta = lattice_fixture.get_reference_particle().get_beta()
    chrom = synergia.simulation.Lattice_simulator.get_chromaticities(lattice_fixture)
    print("H chromaticitity: ", chrom.horizontal_chromaticity)
    print("V chromaticitity: ", chrom.vertical_chromaticity)
    # MADX chromaticiticity is with respect to PT
    # assert chrom.horizontal_chromaticity == pytest.approx(-65.68223681*beta)
    # assert chrom.vertical_chromaticity == pytest.approx(-25.5636675*beta)


def main():
    lf = lattice_fixture()
    test_tunes(lf)
    test_chromaticities(lf)


if __name__ == "__main__":
    main()
