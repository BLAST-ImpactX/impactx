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
        lp[i, 1] = i * 1.0e-4
        lp[i, 3] = (np-i) * 1.0e-4
    bunch.checkin_particles()
    return bunch


def test_adjust_ref_energy1(bunch_fixture):

    refpart = bunch_fixture.get_design_reference_particle()
    energy = refpart.get_total_energy()
    old_p = refpart.get_momentum()

    # add some energy to each particle by changing dp/p.
    dE = 0.050 # 50 MeV
    # what is that in dp/p?
    new_energy = energy + dE
    print('new energy: ', new_energy)
    new_p = np.sqrt(new_energy**2 - mp**2)
    dpop = new_p/old_p - 1.0
    print('dp/p corresponding to new energy: ', dpop)

    # set the dp/p of all the particles
    bunch_fixture.checkout_particles()
    lp = bunch_fixture.get_particles_numpy()
    lp[:, 5] = dpop
    bunch_fixture.checkin_particles()

    # adjust the reference energy
    bunch_fixture.adjust_bunch_particles_reference_energy(new_energy)

    # check that dp/p is set correctly
    bunch_fixture.checkout_particles()
    lp = bunch_fixture.get_particles_numpy()
    print('particles dp/p: ', lp[:, 5])
    for i in range(lp.shape[0]):
        assert lp[i, 5] == pytest.approx(0.0, abs=1.0e-12)
    for i in range(lp.shape[0]):
        assert lp[i, 1] == pytest.approx(i * 1.0e-4 * old_p/new_p)
        assert lp[i, 3] == pytest.approx((lp.shape[0] - i) * 1.0e-4 * old_p/new_p)

def main():
    bf = bunch_fixture()
    test_adjust_ref_energy1(bf)

if __name__ == "__main__":
    main()
