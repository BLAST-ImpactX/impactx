#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 40  # 8 particles in transverse ring x 5 longitudinal positions
realparticles = 4.0e10

kick_values = np.array([32e-4, 16e-4, 8e-4, 4e-4, 2e-4, 1e-6])


# prop_fixture is a propagator
@pytest.fixture
def pf():
    # booster-like lattice
    channel_madx = """
beam, particle=proton,energy=pmass+0.4;

m: matrix, kick1=32e-4, kick2=16e-4, kick3=8e-4, kick4=4e-4, kick5=2e-4, kick6=1e-6;

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


def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        ref_part, macroparticles, realparticles
    )
    bunch = sim.get_bunch()
    s2o2 = np.sqrt(2.0) / 2
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

    dr = 0.001

    # populate longitudinal coordinates first
    k = 0
    for iz in range(-2, 3):
        # loop over longitudinal position
        cdt = 0.5 * iz
        for j in range(8):
            lp[k, 0:6] = 0.0
            lp[k, 4] = cdt
            k = k + 1

    # populate transverse coordinates
    k = 0
    for j in range(5):
        lp[k, 0] = dr
        k = k + 1

        lp[k, 0] = s2o2 * dr
        lp[k, 2] = s2o2 * dr
        k = k + 1

        lp[k, 2] = dr
        k = k + 1

        lp[k, 0] = -s2o2 * dr
        lp[k, 2] = s2o2 * dr
        k = k + 1

        lp[k, 0] = -dr
        k = k + 1

        lp[k, 0] = -s2o2 * dr
        lp[k, 2] = -s2o2 * dr
        k = k + 1

        lp[k, 2] = -dr
        k = k + 1

        lp[k, 0] = s2o2 * dr
        lp[k, 2] = -s2o2 * dr
        k = k + 1

    bunch.checkin_particles()
    return sim


def test_lattice_prop(pf):
    refpart = pf.get_lattice().get_reference_particle()
    sim = create_simulator(refpart)

    bunch = sim.get_bunch(0, 0)
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

    # copy original particles
    op = np.array(lp)[:, :6]

    simlog = synergia.utils.Logger()
    # propagate 1 turn
    pf.propagate(sim, simlog, 1)

    bunch.checkout_particles()
    prop_particles = bunch.get_particles_numpy()

    for i in range(40):
        for j in range(6):
            # print(i, j)
            assert prop_particles[i, j] == pytest.approx(op[i, j] + kick_values[j])


def main():
    pf = pf()
    test_lattice_prop(pf)


if __name__ == "__main__":
    main()
