#!/usr/bin/env python

# this tests that the normal form can be created and applied
# for a simplified Booster model.

import pytest
import numpy
import synergia

simple_booster_src = """
ke1 = 0.8;  !800 MeV kinetic energy at injection

rhof  :=  40.847086;   !  bending radius of focusing magnet
rhod  :=  48.034101;   !  bending radius of defocusing magnet

blength :=     2.889612;    !  arc length for both F and D magnets
blengthf :=    2.889009499; !  physical length (straight length) for F magnet
blengthd :=    2.889176299; !  physical length (straight length) for D magnet

qsd := -57.38855012e-3;
qsf := 54.10921561e-3;

ssd :=  -0.04381647074 + ke1*(0.009150934932+ ke1*(-0.0023900895  + ke1*(0.000318068028 -  ke1* 1.6353205e-05)));
ssf :=  -0.006384940088 + ke1*(0.01967542848 + ke1*( -0.006776746 + ke1*(0.00091367565 - ke1* 4.293705e-05)));

FMAG: SBEND,  L = blength  , ANGLE = blength/rhof, e1=blength/(2*rhof), e2=blength/(2*rhof), K1 = qsf  , K2 = ssf;
DMAG: SBEND,  L = blength  , ANGLE = blength/rhod, e1=blength/(2*rhod), e2=blength/(2*rhod), K1 = qsd  , K2 = ssd;

d1: drift, l=0.6;
d2: drift, l=0.5;
d3: drift, l=6.0;
d4: drift, l=0.5;
d5: drift, l=0.6;

! these drifts pad the RF cavities right after dmagunn beginning the long straight

drlnna: drift, l=0.21;  
drlnnb: drift, l=0.12;
drlnnc: drift, l=0.25;
drlnnd: drift, l=0.111;
cplnn: drift, l=0.168; ! corrector package
drlnne: drift, l=0.251;
hpnnl: drift, l=0.095; ! horizontal bpm
vpnnl: drift, l=0.095; ! vertical bpm

drrf: drift, l=2.35/2;
rfc: rfcavity, l=0, harmon=84, volt=0.2/24 , lag=0;

fmagu14: fmag;
fmagd14: fmag;
dmagu14: dmag;
dmagd14: dmag;
rf01: line=(drrf, rfc, drrf);
rf02: line=(drrf, rfc, drrf);

cell14 : line = (d1, fmagu14, d2, dmagu14, drlnna, rf01, drlnnb, rf02, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd14, d4, fmagd14, d5);

booster: line=(24*cell14);

beam, particle=proton, energy=0.8+pmass;
"""


def test_booster_normal_form():
    reader = synergia.lattice.MadX_reader()
    reader.parse(simple_booster_src)
    lattice = reader.get_lattice("booster")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    # normal form
    # verifying that the normal form can be created.
    # This would fail before PR#137 to fix the trigon log calculation
    # because conjugate eigenvalues could not be found.
    nf = synergia.simulation.Lattice_simulator.calculate_normal_form_o3(lattice)


if __name__ == "__main__":
    test_booster_normal_form()
