#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest


def get_lattice():
    fodo_madx = """
! Simplified Booster lattice

// From JFO 2022-12-08
// The simplified booster lattice is trivial. I have no madx lattice (you could make one
// very easily) - I just use pyorbit classes directly to instantiate a basic cell;  it is then
// replicated 24 times. I took the bending magnets lengths and 
// strengths directly from the official MADX  lattice file. 

// The basic cell is 

// d1 fmag d2 dmag d3 

// d1, d2, d3 : drifts of lengths 0.6 0.5 and 3.0 m
// fmag:  focusing      bend   L = 2.889612 m 
// dmag   defocusing bend   L = 2.889612 m

// total cell length: 19.758 m
// total ring  length = 24*19.758 = 474.20 m 

// The length, focusing strengths and curvature radius of the 
// magnets are as in the booster MADX file.  

// If you entered 1 cell correctly, you should get the periodic solution:
// bx = 33.86 m ax = 0
// by = 5.39m    ay =0 

// For 24 cells, the raw tunes are nux = 7.017 and nuy = 6.675.  You will need to tweak the nominal focusing strengths a bit to avoid resonances.


//--------- Nominal Gradient Magnet Definitions  

// EGS
// The apparent cell structure is actually:

// D1, l=0.6;
// FMAGU01;
// D2, l=0.5;
// DMAGU01;
// D3, l=6.0;
// DMAGD01;
// D4, l=0.5;
// FMAGD01;
// DR, l=0.6



ke1 = 0.8;  !800 MeV kinetic energy at injection

rhof  :=  40.847086;   !  bending radius of focusing magnet
rhod  :=  48.034101;   !  bending radius of defocusing magnet

blength :=     2.889612;    !  arc length for both F and D magnets
blengthf :=    2.889009499; !  physical length (straight length) for F magnet
blengthd :=    2.889176299; !  physical length (straight length) for D magnet


!
! The quad field for the gradient magnet is contained in file " qsdqsf.dat" to be read in before this file !
!
! read from file at time step = 7
qsd := -57.38855012e-3;
qsf := 54.10921561e-3;


! These ssd and ssf strengths come from fitting to 01 Dec 2015 chromaticity data
! and predicts chromaticity much better than using the Drozhdin et al measurements above

ssd :=  -0.04381647074 + ke1*(0.009150934932+ ke1*(-0.0023900895  + ke1*(0.000318068028 -  ke1* 1.6353205e-05)));

ssf :=  -0.006384940088 + ke1*(0.01967542848 + ke1*( -0.006776746 + ke1*(0.00091367565 - ke1* 4.293705e-05)));

 !
 ! Gradient magnets defined by their physical length aith their bend angle
 ! being defined by the arc length/radius of curvature

!FMAG: RBEND,  L = blengthf  , ANGLE = blength/rhof, K1 = qsf  , K2 = ssf;
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
rfc: rfcavity, l=0, harmon=84, volt=0, lag=0;

!!!!!!!!!!!!!!!!   beginning of ring definition

fmagu01: fmag;
fmagd01: fmag;
dmagu01: dmag;
dmagd01: dmag;

cell01 : line = (d1, fmagu01, d2, dmagu01, d3, dmagd01, d4, fmagd01, d5);

fmagu02: fmag;
fmagd02: fmag;
dmagu02: dmag;
dmagd02: dmag;

cell02 : line = (d1, fmagu02, d2, dmagu02, d3, dmagd02, d4, fmagd02, d5);

fmagu03: fmag;
fmagd03: fmag;
dmagu03: dmag;
dmagd03: dmag;

cell03 : line = (d1, fmagu03, d2, dmagu03, d3, dmagd03, d4, fmagd03, d5);

fmagu04: fmag;
fmagd04: fmag;
dmagu04: dmag;
dmagd04: dmag;

cell04 : line = (d1, fmagu04, d2, dmagu04, d3, dmagd04, d4, fmagd04, d5);

fmagu05: fmag;
fmagd05: fmag;
dmagu05: dmag;
dmagd05: dmag;

cell05 : line = (d1, fmagu05, d2, dmagu05, d3, dmagd05, d4, fmagd05, d5);

fmagu06: fmag;
fmagd06: fmag;
dmagu06: dmag;
dmagd06: dmag;

cell06 : line = (d1, fmagu06, d2, dmagu06, d3, dmagd06, d4, fmagd06, d5);

fmagu07: fmag;
fmagd07: fmag;
dmagu07: dmag;
dmagd07: dmag;

cell07 : line = (d1, fmagu07, d2, dmagu07, d3, dmagd07, d4, fmagd07, d5);

fmagu08: fmag;
fmagd08: fmag;
dmagu08: dmag;
dmagd08: dmag;

cell08 : line = (d1, fmagu08, d2, dmagu08, d3, dmagd08, d4, fmagd08, d5);

fmagu09: fmag;
fmagd09: fmag;
dmagu09: dmag;
dmagd09: dmag;

cell09 : line = (d1, fmagu09, d2, dmagu09, d3, dmagd09, d4, fmagd09 , d5);

fmagu10: fmag;
fmagd10: fmag;
dmagu10: dmag;
dmagd10: dmag;

cell10 : line = (d1, fmagu10, d2, dmagu10, d3, dmagd10, d4, fmagd10 , d5);

fmagu11: fmag;
fmagd11: fmag;
dmagu11: dmag;
dmagd11: dmag;

cell11 : line = (d1, fmagu11, d2, dmagu11, d3, dmagd11, d4, fmagd11, d5);

fmagu12: fmag;
fmagd12: fmag;
dmagu12: dmag;
dmagd12: dmag;

cell12 : line = (d1, fmagu12, d2, dmagu12, d3, dmagd12, d4, fmagd12, d5);

fmagu13: fmag;
fmagd13: fmag;
dmagu13: dmag;
dmagd13: dmag;

cell13 : line = (d1, fmagu13, d2, dmagu13, d3, dmagd13, d4, fmagd13, d5);

fmagu14: fmag;
fmagd14: fmag;
dmagu14: dmag;
dmagd14: dmag;
rf01: line=(drrf, rfc, drrf);
rf02: line=(drrf, rfc, drrf);

cell14 : line = (d1, fmagu14, d2, dmagu14, drlnna, rf01, drlnnb, rf02, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd14, d4, fmagd14, d5);

fmagu15: fmag;
fmagd15: fmag;
dmagu15: dmag;
dmagd15: dmag;
rf03: line=(drrf, rfc, drrf);
rf04: line=(drrf, rfc, drrf);

cell15 : line = (d1, fmagu15, d2, dmagu15, drlnna, rf03, drlnnb, rf04, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd15, d4, fmagd15, d5);

fmagu16: fmag;
fmagd16: fmag;
dmagu16: dmag;
dmagd16: dmag;
rf05: line=(drrf, rfc, drrf);
rf06: line=(drrf, rfc, drrf);

cell16 : line = (d1, fmagu16, d2, dmagu16, drlnna, rf05, drlnnb, rf06, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd16, d4, fmagd16, d5);

fmagu17: fmag;
fmagd17: fmag;
dmagu17: dmag;
dmagd17: dmag;
rf07: line=(drrf, rfc, drrf);
rf08: line=(drrf, rfc, drrf);

cell17 : line = (d1, fmagu17, d2, dmagu17, drlnna, rf07, drlnnb, rf08, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd17, d4, fmagd17, d5);

fmagu18: fmag;
fmagd18: fmag;
dmagu18: dmag;
dmagd18: dmag;
rf09: line=(drrf, rfc, drrf);
rf10: line=(drrf, rfc, drrf);

cell18 : line = (d1, fmagu18, d2, dmagu18, drlnna, rf09, drlnnb, rf10, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd18, d4, fmagd18, d5);

fmagu19: fmag;
fmagd19: fmag;
dmagu19: dmag;
dmagd19: dmag;
rf11: line=(drrf, rfc, drrf);
rf12: line=(drrf, rfc, drrf);

cell19 : line = (d1, fmagu19, d2, dmagu19, drlnna, rf11, drlnnb, rf12, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd19, d4, fmagd19 , d5);

fmagu20: fmag;
fmagd20: fmag;
dmagu20: dmag;
dmagd20: dmag;
rf13: line=(drrf, rfc, drrf);
rf14: line=(drrf, rfc, drrf);

cell20 : line = (d1, fmagu20, d2, dmagu20, drlnna, rf13, drlnnb, rf14, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd20, d4, fmagd20 , d5);

fmagu21: fmag;
fmagd21: fmag;
dmagu21: dmag;
dmagd21: dmag;
rf15: line=(drrf, rfc, drrf);
rf16: line=(drrf, rfc, drrf);

cell21 : line = (d1, fmagu21, d2, dmagu21, drlnna, rf15, drlnnb, rf16, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd21, d4, fmagd21, d5);

fmagu22: fmag;
fmagd22: fmag;
dmagu22: dmag;
dmagd22: dmag;
rf17: line=(drrf, rfc, drrf);
rf18: line=(drrf, rfc, drrf);

cell22 : line = (d1, fmagu22, d2, dmagu22, drlnna, rf17, drlnnb, rf18, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd22, d4, fmagd22, d5);

fmagu23: fmag;
fmagd23: fmag;
dmagu23: dmag;
dmagd23: dmag;
rf19: line=(drrf, rfc, drrf);
rf20: line=(drrf, rfc, drrf);

cell23 : line = (d1, fmagu23, d2, dmagu23, drlnna, rf19, drlnnb, rf20, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd23, d4, fmagd23, d5);

fmagu24: fmag;
fmagd24: fmag;
dmagu24: dmag;
dmagd24: dmag;
rf21: line=(drrf, rfc, drrf);
rf22: line=(drrf, rfc, drrf);

cell24 : line = (d1, fmagu24, d2, dmagu24, drlnna, rf21, drlnnb, rf22, drlnnc, hpnnl, vpnnl, drlnnd, cplnn, drlnne, dmagd24, d4, fmagd24, d5);

booster: line=(cell01, cell02, cell03, cell04, cell05, cell06, cell07, cell08, cell09, cell10, cell11, cell12, cell13, cell14, cell15, cell16, cell17, cell18, cell19, cell20, cell21, cell22, cell23, cell24);

!beam, particle=proton, energy=0.8+pmass;
beam, particle=proton, energy=1.7388477415186678; ! this energy causes closed orbit failures at a tolerance of 1.0e-13

"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(fodo_madx)
    lattice = reader.get_lattice("booster")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    return lattice


def test_correct_total_energy():
    lattice = get_lattice()
    assert lattice.get_reference_particle().get_total_energy() == pytest.approx(
        1.7388477415186678
    )
    assert lattice.get_lattice_energy() == pytest.approx(1.7388477415186678)


def test_correct_length():
    lattice = get_lattice()
    assert lattice.get_length() == pytest.approx(474.202752)


def test_correct_init_frequency():
    lattice = get_lattice()
    length = lattice.get_length()
    # harmonic number is unity
    freq = (
        84
        * lattice.get_reference_particle().get_beta()
        * synergia.foundation.pconstants.c
        / length
    )
    ncavities = 0
    rfelem = None
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            ncavities = ncavities + 1
            rfelem = elem
    assert ncavities == 22
    rffreq = rfelem.get_double_attribute("freq")
    # rffreq will be in MHz, while my calculated freq is in Hz
    assert rffreq * 1.0e6 == pytest.approx(freq)


def test_correct_cdt():
    lattice = get_lattice()
    beta = lattice.get_reference_particle().get_beta()
    length = lattice.get_length()
    # get closed orbit state from tune_circular_lattice
    closed_orbit_state = synergia.simulation.Lattice_simulator.tune_circular_lattice(
        lattice
    )
    cdt = closed_orbit_state[4]
    assert length / beta == pytest.approx(cdt)


def test_correct_increased_frequency():
    lattice = get_lattice()
    length = lattice.get_length()
    # harmonic number is unity
    init_energy = lattice.get_reference_particle().get_total_energy()
    print("init_energy: ", init_energy)
    gamma = init_energy / synergia.foundation.pconstants.mp
    betagamma = np.sqrt(gamma**2 - 1)
    beta = betagamma / gamma
    print("orig betagamma: ", betagamma)
    print("orig gamma: ", gamma)
    print("orig beta: ", beta)
    # return
    new_energy = init_energy + 1
    print("new_energy: ", new_energy)
    # lattice.get_reference_particle().set_total_energy(new_energy)
    lattice.set_lattice_energy(new_energy)
    # return
    print("readback new energy: ", lattice.get_reference_particle().get_total_energy())
    new_gamma = new_energy / synergia.foundation.pconstants.mp
    new_betagamma = np.sqrt(new_gamma**2 - 1)
    new_beta = new_betagamma / new_gamma
    print("new_betagamma: ", new_betagamma)
    print("new_gamma: ", new_gamma)
    print("readback new gamma: ", lattice.get_reference_particle().get_gamma())
    print("new_beta: ", new_beta)
    print("readback new beta: ", lattice.get_reference_particle().get_beta())
    assert lattice.get_reference_particle().get_beta() == pytest.approx(new_beta)
    print("before tune_circular_lattice", flush=True)
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    print("after tune_circular_lattice", flush=True)
    print(
        "readback beta after tune_circular_lattice: ",
        lattice.get_reference_particle().get_beta(),
    )
    print("readback energy after tune_circular_lattice: ", lattice.get_lattice_energy())
    freq = (
        84
        * lattice.get_reference_particle().get_beta()
        * synergia.foundation.pconstants.c
        / length
    )
    ncavities = 0
    rfelem = None
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            ncavities = ncavities + 1
            rfelem = elem
    assert ncavities == 22
    rffreq = rfelem.get_double_attribute("freq")
    # rffreq will be in MHz, while my calculated freq is in Hz
    print("rffreq: ", rffreq * 1.0e6)
    print("freq: ", freq)
    assert rffreq * 1.0e6 == pytest.approx(freq)


def main():
    print(get_lattice())
    test_correct_length()
    test_correct_total_energy()
    test_correct_init_frequency()
    test_correct_increased_frequency()


if __name__ == "__main__":
    main()
