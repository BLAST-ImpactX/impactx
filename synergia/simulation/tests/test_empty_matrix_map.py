#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10

nturns = 100


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    # booster-like lattice
    channel_madx = """
beam, particle=proton,energy=pmass+0.8;
! totally empty matrix element should have identity map
m: matrix;
channel: sequence, l=0.0;
m, at=0.0;
endsequence;

"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
    print(lattice)
    lattice.set_all_string_attribute("extractor_type", "libff")
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


def test_lattice_map(prop_fixture):
    lattice = prop_fixture.get_lattice()
    map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)

    # diagonal elements should be 1, all others 0
    for i in range(6):
        for j in range(6):
            if i == j:
                assert map[i, j] == pytest.approx(1.0)
            else:
                assert map[i, j] == pytest.approx(0.0)


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


def main():
    pf = prop_fixture()
    test_lattice_map(pf)


if __name__ == "__main__":
    main()
