#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

import synergia.simulation as SIM
ET = synergia.lattice.element_type
MT = synergia.lattice.marker_type
PCONST = synergia.foundation.pconstants

kinetic_energy = 0.00250 # IOTA Kinetic Energy 2.5 MeV
stepper = "elements"
steps = 1
harmonic_number = 4

################################################################################

iota_sequence = """
kqa1r := kq01;
kq01 = -7.72652301;
kqa2r := kq02;
kq02 = 12.28401222;
kqa3r := kq03;
kq03 = -12.43016989;
kqa4r := kq04;
kq04 = 20.16074347;
kqb1r := kq05;
kq05 = -10.24365752;
kqb2r := kq06;
kq06 = 15.12808788;
kqb3r := kq07;
kq07 = -6.92311681;
kqb4r := kq08;
kq08 = -6.90057605;
kqb5r := kq09;
kq09 = 13.50655178;
kqb6r := kq10;
kq10 = -11.91343344;
kqc1r := kq11;
kq11 = -13.51948869;
kqc2r := kq12;
kq12 = 12.0339278;
kqc3r := kq13;
kq13 = -13.56878135;
kqd1r := kq14;
kq14 = -7.97007816;
kqd2r := kq15;
kq15 = 5.92322639;
kqd3r := kq16;
kq16 = -6.32915747;
kqd4r := kq17;
kq17 = 5.1636516;
kqe1r := kq18;
kq18 = -4.69712477;
kqe2r := kq19;
kq19 = 7.0326898;
kqe3 := kq20;
kq20 = -7.19881671;
kqe2l := kq19;
kqe1l := kq18;
kqd4l := kq17;
kqd3l := kq16;
kqd2l := kq15;
kqd1l := kq14;
kqc3l := kq13;
kqc2l := kq12;
kqc1l := kq11;
kqb6l := kq10;
kqb5l := kq09;
kqb4l := kq08;
kqb3l := kq07;
kqb2l := kq06;
kqb1l := kq05;
kqa4l := kq04;
kqa3l := kq03;
kqa2l := kq02;
kqa1l := kq01;
mseqari: marker;
ibpm: monitor;
ibpma1: ibpm;
qa1r: quadrupole,l:= 0.21,k1:=kqa1r ;
qa2r: quadrupole,l:= 0.21,k1:=kqa2r ;
sqa1r: quadrupole,l:= 0.1,k1s:= 0;
ibpma2r: ibpm;
qa3r: quadrupole,l:= 0.21,k1:=kqa3r ;
qa4r: quadrupole,l:= 0.21,k1:=kqa4r ;
ibpma3r: ibpm;
sqa2r: quadrupole,l:= 0.1,k1s:= 0;
mseqare: marker;
mphm1ri: marker;
dedge30: dipedge,e1:= 0,h:= 1.338646717,hgap:= 0.015042,fint:= 0.5;
m1r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm1re: marker;
mseqbri: marker;
sqb1r: quadrupole,l:= 0.1,k1s:= 0;
qb1r: quadrupole,l:= 0.21,k1:=kqb1r ;
qb2r: quadrupole,l:= 0.21,k1:=kqb2r ;
qb3r: quadrupole,l:= 0.21,k1:=kqb3r ;
ibpmb1r: ibpm;
nlr1: marker;
ior: marker;
nlr2: marker;
ibpmb2r: ibpm;
qb4r: quadrupole,l:= 0.21,k1:=kqb4r ;
qb5r: quadrupole,l:= 0.21,k1:=kqb5r ;
qb6r: quadrupole,l:= 0.21,k1:=kqb6r ;
sqb2r: quadrupole,l:= 0.1,k1s:= 0;
mseqbre: marker;
mphm2ri: marker;
dedge60: dipedge,e1:= 0,h:= 1.381554029,hgap:= 0.014786,fint:= 0.5;
m2r: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm2re: marker;
mseqcri: marker;
sqc1r: quadrupole,l:= 0.1,k1s:= 0;
ibpmc1r: ibpm;
qc1r: quadrupole,l:= 0.21,k1:=kqc1r ;
qc2r: quadrupole,l:= 0.21,k1:=kqc2r ;
qc3r: quadrupole,l:= 0.21,k1:=kqc3r ;
ibpmc2r: ibpm;
sqc2r: quadrupole,l:= 0.1,k1s:= 0;
mseqcre: marker;
mphm3ri: marker;
m3r: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm3re: marker;
mseqdri: marker;
ibpmd1r: ibpm;
sqd1r: quadrupole,l:= 0.1,k1s:= 0;
qd1r: quadrupole,l:= 0.21,k1:=kqd1r ;
qd2r: quadrupole,l:= 0.21,k1:=kqd2r ;
el1: marker;
cel: solenoid,l:= 0.7,ks:= 0;
el2: marker;
qd3r: quadrupole,l:= 0.21,k1:=kqd3r ;
sqd2r: quadrupole,l:= 0.1,k1s:= 0;
qd4r: quadrupole,l:= 0.21,k1:=kqd4r ;
ibpmd2r: ibpm;
mseqdre: marker;
mphm4ri: marker;
m4r: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm4re: marker;
mseqei: marker;
ibpme1r: ibpm;
qe1r: quadrupole,l:= 0.21,k1:=kqe1r ;
sqe1r: quadrupole,l:= 0.1,k1s:= 0;
ibpme2r: ibpm;
qe2r: quadrupole,l:= 0.21,k1:=kqe2r ;
sqe2r: quadrupole,l:= 0.1,k1s:= 0;
qe3: quadrupole,l:= 0.21,k1:=kqe3 ;
sqe2l: quadrupole,l:= 0.1,k1s:= 0;
qe2l: quadrupole,l:= 0.21,k1:=kqe2l ;
ibpme2l: ibpm;
sqe1l: quadrupole,l:= 0.1,k1s:= 0;
qe1l: quadrupole,l:= 0.21,k1:=kqe1l ;
ibpme1l: ibpm;
mseqee: marker;
mphm4li: marker;
m4l: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm4le: marker;
mseqdli: marker;
ibpmd2l: ibpm;
qd4l: quadrupole,l:= 0.21,k1:=kqd4l ;
sqd2l: quadrupole,l:= 0.1,k1s:= 0;
qd3l: quadrupole,l:= 0.21,k1:=kqd3l ;
rfc: rfcavity,l:= 0.05,volt:= 0.000847,lag:= 0,harmon:= 4;
qd2l: quadrupole,l:= 0.21,k1:=kqd2l ;
qd1l: quadrupole,l:= 0.21,k1:=kqd1l ;
sqd1l: quadrupole,l:= 0.1,k1s:= 0;
ibpmd1l: ibpm;
mseqdle: marker;
mphm3li: marker;
m3l: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm3le: marker;
mseqcli: marker;
sqc2l: quadrupole,l:= 0.1,k1s:= 0;
ibpmc2l: ibpm;
qc3l: quadrupole,l:= 0.21,k1:=kqc3l ;
qc2l: quadrupole,l:= 0.21,k1:=kqc2l ;
qc1l: quadrupole,l:= 0.21,k1:=kqc1l ;
ibpmc1l: ibpm;
sqc1l: quadrupole,l:= 0.1,k1s:= 0;
mseqcle: marker;
mphm2li: marker;
m2l: sbend,l:= 0.757985232,angle:= 1.047197551;
mphm2le: marker;
mseqbli: marker;
sqb2l: quadrupole,l:= 0.1,k1s:= 0;
qb6l: quadrupole,l:= 0.21,k1:=kqb6l ;
qb5l: quadrupole,l:= 0.21,k1:=kqb5l ;
qb4l: quadrupole,l:= 0.21,k1:=kqb4l ;
ibpmb2l: ibpm;
nll1: marker;
nll2: marker;
ibpmb1l: ibpm;
qb3l: quadrupole,l:= 0.21,k1:=kqb3l ;
qb2l: quadrupole,l:= 0.21,k1:=kqb2l ;
qb1l: quadrupole,l:= 0.21,k1:=kqb1l ;
sqb1l: quadrupole,l:= 0.1,k1s:= 0;
mseqble: marker;
mphm1li: marker;
m1l: sbend,l:= 0.3911403725,angle:= 0.5235987756;
mphm1le: marker;
mseqali: marker;
sqa2l: quadrupole,l:= 0.1,k1s:= 0;
ibpma3l: ibpm;
qa4l: quadrupole,l:= 0.21,k1:=kqa4l ;
qa3l: quadrupole,l:= 0.21,k1:=kqa3l ;
ibpma2l: ibpm;
sqa1l: quadrupole,l:= 0.1,k1s:= 0;
qa2l: quadrupole,l:= 0.21,k1:=kqa2l ;
qa1l: quadrupole,l:= 0.21,k1:=kqa1l ;
mseqale: marker;
iota: sequence, l = 39.95567226;
mseqari, at = 0;
ibpma1, at = 0.02;
qa1r, at = 1.0175;
qa2r, at = 1.3625;
sqa1r, at = 1.99;
ibpma2r, at = 2.095;
qa3r, at = 2.2975;
qa4r, at = 2.6525;
ibpma3r, at = 2.865;
sqa2r, at = 2.97;
mseqare, at = 3.0405;
mphm1ri, at = 3.0405;
dedge30, at = 3.117400202;
m1r, at = 3.312970389;
dedge30, at = 3.508540575;
mphm1re, at = 3.585440777;
mseqbri, at = 3.585440777;
sqb1r, at = 3.715940777;
qb1r, at = 3.953440777;
qb2r, at = 4.303440777;
qb3r, at = 4.653440777;
ibpmb1r, at = 4.865940777;
nlr1, at = 4.910940777;
ior, at = 5.810940777;
nlr2, at = 6.710940777;
ibpmb2r, at = 6.755940777;
qb4r, at = 6.968440777;
qb5r, at = 7.318440777;
qb6r, at = 7.668440777;
sqb2r, at = 7.905940777;
mseqbre, at = 8.036440777;
mphm2ri, at = 8.036440777;
dedge60, at = 8.112186805;
m2r, at = 8.491179421;
dedge60, at = 8.870172037;
mphm2re, at = 8.945918065;
mseqcri, at = 8.945918065;
sqc1r, at = 9.106418065;
ibpmc1r, at = 9.211418065;
qc1r, at = 9.423918065;
qc2r, at = 9.988918065;
qc3r, at = 10.55391806;
ibpmc2r, at = 10.76641806;
sqc2r, at = 10.87141806;
mseqcre, at = 11.03191806;
mphm3ri, at = 11.03191806;
dedge60, at = 11.10766409;
m3r, at = 11.48665671;
dedge60, at = 11.86564932;
mphm3re, at = 11.94139535;
mseqdri, at = 11.94139535;
ibpmd1r, at = 12.20689535;
sqd1r, at = 12.39189535;
qd1r, at = 12.61939535;
qd2r, at = 13.24939535;
el1, at = 13.81689535;
cel, at = 14.16689535;
el2, at = 14.51689535;
qd3r, at = 15.08439535;
sqd2r, at = 15.39939535;
qd4r, at = 15.71439535;
ibpmd2r, at = 16.12689535;
mseqdre, at = 16.39239535;
mphm4ri, at = 16.39239535;
dedge30, at = 16.46929555;
m4r, at = 16.66486574;
dedge30, at = 16.86043593;
mphm4re, at = 16.93733613;
mseqei, at = 16.93733613;
ibpme1r, at = 17.20283613;
qe1r, at = 17.51533613;
sqe1r, at = 17.76283613;
ibpme2r, at = 18.74783613;
qe2r, at = 18.96033613;
sqe2r, at = 19.34783613;
qe3, at = 19.97783613;
sqe2l, at = 20.60783613;
qe2l, at = 20.99533613;
ibpme2l, at = 21.20783613;
sqe1l, at = 22.19283613;
qe1l, at = 22.44033613;
ibpme1l, at = 22.75283613;
mseqee, at = 23.01833613;
mphm4li, at = 23.01833613;
dedge30, at = 23.09523633;
m4l, at = 23.29080652;
dedge30, at = 23.4863767;
mphm4le, at = 23.56327691;
mseqdli, at = 23.56327691;
ibpmd2l, at = 23.82877691;
qd4l, at = 24.24127691;
sqd2l, at = 24.55627691;
qd3l, at = 24.87127691;
rfc, at = 25.78877691;
qd2l, at = 26.70627691;
qd1l, at = 27.33627691;
sqd1l, at = 27.56377691;
ibpmd1l, at = 27.74877691;
mseqdle, at = 28.01427691;
mphm3li, at = 28.01427691;
dedge60, at = 28.09002293;
m3l, at = 28.46901555;
dedge60, at = 28.84800817;
mphm3le, at = 28.92375419;
mseqcli, at = 28.92375419;
sqc2l, at = 29.08425419;
ibpmc2l, at = 29.18925419;
qc3l, at = 29.40175419;
qc2l, at = 29.96675419;
qc1l, at = 30.53175419;
ibpmc1l, at = 30.74425419;
sqc1l, at = 30.84925419;
mseqcle, at = 31.00975419;
mphm2li, at = 31.00975419;
dedge60, at = 31.08550022;
m2l, at = 31.46449284;
dedge60, at = 31.84348545;
mphm2le, at = 31.91923148;
mseqbli, at = 31.91923148;
sqb2l, at = 32.04973148;
qb6l, at = 32.28723148;
qb5l, at = 32.63723148;
qb4l, at = 32.98723148;
ibpmb2l, at = 33.19973148;
nll1, at = 33.24473148;
nll2, at = 35.04473148;
ibpmb1l, at = 35.08973148;
qb3l, at = 35.30223148;
qb2l, at = 35.65223148;
qb1l, at = 36.00223148;
sqb1l, at = 36.23973148;
mseqble, at = 36.37023148;
mphm1li, at = 36.37023148;
dedge30, at = 36.44713168;
m1l, at = 36.64270187;
dedge30, at = 36.83827206;
mphm1le, at = 36.91517226;
mseqali, at = 36.91517226;
sqa2l, at = 36.98567226;
ibpma3l, at = 37.09067226;
qa4l, at = 37.30317226;
qa3l, at = 37.65817226;
ibpma2l, at = 37.86067226;
sqa1l, at = 37.96567226;
qa2l, at = 38.59317226;
qa1l, at = 38.93817226;
mseqale, at = 39.95567226;
endsequence;
beam, particle=proton, energy=pmass+0.00250;
"""


################################################################################


@pytest.fixture
def lattice_fixture():
    # read the lattice in from a MadX sequence file

    reader = synergia.lattice.MadX_reader()
    reader.parse(iota_sequence)
    lattice = reader.get_lattice("iota")

    # The sequence doesn't have a reference particle so define it here
    KE = kinetic_energy
    mass = PCONST.mp

    etot = KE + mass
    refpart = synergia.foundation.Reference_particle(1, mass, etot)

    lattice.set_reference_particle(refpart)

    return lattice

################################################################################

def get_propagator(lattice):

    sc_ops = synergia.collective.Dummy_CO_options()
    stepper = synergia.simulation.Independent_stepper_elements(steps)
    
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

################################################################################

def test_elem_lengths(lattice_fixture):
    lattice = lattice_fixture
    lattice_length = lattice.get_length()

    refpart = lattice.get_reference_particle()
    beta = refpart.get_beta()
    
    # set the RF cavities
    state = SIM.Lattice_simulator.tune_circular_lattice(lattice)
    
    # the cdt of state should be the orbit length/beta
    assert lattice_length == pytest.approx(state[4]*beta)

    # Check the frequency of the RF cavities
    freq_should_be = harmonic_number * beta * PCONST.c/lattice_length

    for elem in lattice.get_elements():
        # is this an RF cavity?
        if elem.get_type() == ET.rfcavity:
            # RF cavity freq is in MHz
            freq = elem.get_double_attribute('freq')*1.0e6
            assert freq_should_be == pytest.approx(freq)

    # Propagate a bunch and check each element length
    element_lengths = [e.get_length() for e in lattice.get_elements()]

    # logger for simulation
    simlog = synergia.utils.parallel_utils.Logger(0, 
            synergia.utils.parallel_utils.LoggerV.INFO_TURN)

    commxx = synergia.utils.Commxx()
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        refpart, 16, 1.0e10, commxx)

    bunch = sim.get_bunch(0, 0)
    parts = bunch.get_particles_numpy()
    parts[:, 0:6] = 0.0
    bunch.checkin_particles()

    # use step_end_action to capture the cdt for each element
    class context:
        elem_cdt = []

    def step_end_action(sim, inlattice, turn, step):
        bunch = sim.get_bunch(0, 0)
        cdt = bunch.get_design_reference_particle().get_state()[4]
        context.elem_cdt.append(cdt)
        return
    
    sim.reg_prop_action_step_end(step_end_action)

    propagator = get_propagator(lattice)
    propagator.propagate(sim, simlog, 1)

    assert len(context.elem_cdt) == len(element_lengths)

    # Now the elem_cdts * beta should equal the element lengths
    for i in range(len(element_lengths)):
        #print(i, element_lengths[i], context.elem_cdt[i]*beta, flush=True)
        assert element_lengths[i] == pytest.approx(context.elem_cdt[i]*beta)

def main():
    lf = lattice_fixture()
    test_elem_lengths(lf)

if __name__ == "__main__":
    main()
