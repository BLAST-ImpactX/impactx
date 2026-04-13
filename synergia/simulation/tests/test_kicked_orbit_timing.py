#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10
nturns = 10


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    # booster-like lattice
    booster_madx = """
beam, particle=proton,energy=pmass+0.8;
k1: vkicker, vkick=0.01;
k2: vkicker, vkick=-0.02;
k3: vkicker, vkick=0.01;
rfc: rfcavity, l=0.0, volt=0.2, harmon=1;
channel: sequence, l=15, refer=entry;
k1, at=0.0;
k2, at=5.0;
k3, at=10.0;
rfc, at=12.0;
endsequence;

"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(booster_madx)
    lattice = reader.get_lattice("channel")
    lattice.set_all_string_attribute("extractor_type", "libff")
    # The frequency should be decreased to account for the orbit distortion
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


def test_lattice_energy(prop_fixture):
    energy = prop_fixture.get_lattice().get_lattice_energy()
    assert energy == pytest.approx(synergia.foundation.pconstants.mp + 0.8)


def test_lattice_length(prop_fixture):
    assert prop_fixture.get_lattice().get_length() == 15.0


def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        ref_part, macroparticles, realparticles
    )
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch.checkin_particles()
    return sim


def test_proper_frequency(prop_fixture):
    lattice = prop_fixture.get_lattice()
    refpart = lattice.get_reference_particle()
    beta = refpart.get_beta()
    pathlength = 2 * np.sqrt(5**2 + 0.05**2) + 5
    assert pathlength != lattice.get_length()

    # find the rf cavity
    rfc = None
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.rfcavity:
            rfc = elem
            break
    if rfc is None:
        raise RuntimeError("no cavity found")
    cav_freq = rfc.get_double_attribute("freq") * 1.0e6

    print("0 length kick calculated pathlength: ", pathlength)
    print(
        "freq calculated pathlength: ",
        beta * synergia.foundation.pconstants.c / cav_freq,
    )

    assert (
        beta * synergia.foundation.pconstants.c / lattice.get_length()
        != pytest.approx(cav_freq)
    )

    assert beta * synergia.foundation.pconstants.c / pathlength == pytest.approx(
        cav_freq
    )

    # propagate bunch around and make sure it doesn't wobble longitudinally
    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    prop_fixture.propagate(sim, simlog, nturns)

    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

    assert lp[0, 4] == pytest.approx(0.0)
    assert lp[0, 5] == pytest.approx(0.0)


def main():
    pf = prop_fixture()
    test_proper_frequency(pf)


if __name__ == "__main__":
    main()
