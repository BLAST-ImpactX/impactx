#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest


def get_lattice():
    fodo_madx = """
beam, particle=proton,pc=0.75*pmass;

f: quadrupole, l=1.0, k1=0.0625;
d: quadrupole, l=1.0, k1=-0.0625;
rfc: rfcavity, l=0.0, volt=0.2, harmon=1;

fodo: sequence, l=20.0, refer=centre;
fodo_1: f, at=1.0;
fodo_2: d, at=9.0;
fodo_3: d, at=11.0;
fodo_4: f, at=19.0;
fodo_5: rfc, at=20.0;
endsequence;

! beta_x == 32.1571909
! beta_y == 10.3612857
! alphax == alphay == 0
! q_x == q_y == 0.18409
"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(fodo_madx)
    lattice = reader.get_lattice("fodo")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    return lattice


def test_correct_total_energy():
    lattice = get_lattice()
    assert lattice.get_reference_particle().get_total_energy() == pytest.approx(
        1.25 * synergia.foundation.pconstants.mp
    )
    assert lattice.get_lattice_energy() == pytest.approx(
        1.25 * synergia.foundation.pconstants.mp
    )


def test_correct_length():
    lattice = get_lattice()
    assert lattice.get_length() == pytest.approx(20.0)


def test_correct_init_frequency():
    lattice = get_lattice()
    length = lattice.get_length()
    # harmonic number is unity
    freq = (
        lattice.get_reference_particle().get_beta()
        * synergia.foundation.pconstants.c
        / length
    )
    ncavities = 0
    rfelem = None
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            ncavities = ncavities + 1
            rfelem = elem
    assert ncavities == 1
    rffreq = rfelem.get_double_attribute("freq")
    # rffreq will be in MHz, while my calculated freq is in Hz
    assert rffreq * 1.0e6 == pytest.approx(freq)


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
        lattice.get_reference_particle().get_beta()
        * synergia.foundation.pconstants.c
        / length
    )
    ncavities = 0
    rfelem = None
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            ncavities = ncavities + 1
            rfelem = elem
    assert ncavities == 1
    rffreq = rfelem.get_double_attribute("freq")
    # rffreq will be in MHz, while my calculated freq is in Hz
    print("rffreq: ", rffreq * 1.0e6)
    print("freq: ", freq)
    assert rffreq * 1.0e6 == pytest.approx(freq)


def test_correct_increased_frequency2():
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
    return
    new_refpart = lattice.get_reference_particle()
    new_energy = init_energy + 1
    print("new_energy: ", new_energy)
    new_refpart.set_total_energy(new_energy)
    lattice.set_reference_particle(new_refpart)
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
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    freq = (
        lattice.get_reference_particle().get_beta()
        * synergia.foundation.pconstants.c
        / length
    )
    ncavities = 0
    rfelem = None
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            ncavities = ncavities + 1
            rfelem = elem
    assert ncavities == 1
    rffreq = rfelem.get_double_attribute("freq")
    # rffreq will be in MHz, while my calculated freq is in Hz
    assert rffreq * 1.0e6 == pytest.approx(freq)


def main():
    print(get_lattice())
    test_correct_length()
    test_correct_total_energy()
    test_correct_init_frequency()
    test_correct_increased_frequency()


if __name__ == "__main__":
    main()
