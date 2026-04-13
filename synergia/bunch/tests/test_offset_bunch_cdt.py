#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles=16
realparticles=4.0e10

KE0 = 0.8 # starting kinetic energy
mp = synergia.foundation.pconstants.mp

# prop_fixture is a Bunch
@pytest.fixture
def bunch_fixture():
    ref_part = synergia.foundation.Reference_particle(1, mp, mp+KE0)
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(ref_part, macroparticles, realparticles)
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    np = lp.shape[0]
    for i in range(np):
        lp[i, 4] = i * 0.05
    bunch.checkin_particles()
    return bunch


def test_offset_bunch_cdt(bunch_fixture):

    bunch_fixture.checkout_particles()
    lp = bunch_fixture.get_particles_numpy()
    print('particles cdt before offset: ', lp[:, 4])

    # offset the cdt
    bunch_fixture.offset_bunch_particles_cdt(-0.02)

    # check that cdt is set correctly
    # particles changed on device so check them out again
    bunch_fixture.checkout_particles()
    lp = bunch_fixture.get_particles_numpy()
    print('particles cdt after offset: ', lp[:, 4])
    for i in range(lp.shape[0]):
        assert lp[i, 4] == pytest.approx(i*0.05-0.02, abs=1.0e-12)

def main():
    bf = bunch_fixture()
    test_offset_bunch_cdt(bf)

if __name__ == "__main__":
    main()
