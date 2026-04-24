#!/usr/bin/env python
import numpy as np
import pytest

import synergia

macroparticles = 25
realparticles = 4.0e10

nturns = 1


mp = synergia.foundation.pconstants.mp
xoffs = 0.0005
# xpoffs = 0.0001 # need an xp offset to generate edge kick
xpoffs = 0.0
yoffs = 0.0005

##########################################################################################


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture_edgekick():
    channel_madx = """
beam, particle=proton,energy=0.8+pmass;
b: sbend, angle=pi/30.0, l=3.0, e1=pi/60, e2=pi/60;
channel: sequence, l=3.0, refer=entry;
b, at=0.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
    lattice.set_all_string_attribute("extractor_type", "libff")
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


##########################################################################################


# @pytest.fixture
def prop_fixture_noentrykick():
    channel_madx = """
beam, particle=proton,energy=0.8+pmass;
b: sbend, angle=pi/30.0, l=3.0, e1=pi/60.0, e2=pi/60.0, kill_entry_kick=1.0;
channel: sequence, l=3.0, refer=entry;
b, at=0.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
    lattice.set_all_string_attribute("extractor_type", "libff")
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


##########################################################################################


# @pytest.fixture
def prop_fixture_noexitkick():
    channel_madx = """
beam, particle=proton,energy=0.8+pmass;
b: sbend, angle=pi/30.0, l=3.0, e1=pi/60.0, e2=pi/60.0, kill_exit_kick=1.0;
channel: sequence, l=3.0, refer=entry;
b, at=0.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
    lattice.set_all_string_attribute("extractor_type", "libff")
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


##########################################################################################


@pytest.fixture
def prop_fixture_nokick():
    channel_madx = """
beam, particle=proton,energy=0.8+pmass;
b: sbend, angle=pi/30.0, l=3.0, e1=pi/60.0, e2=pi/60.0, kill_entry_kick=1.0, kill_exit_kick=1.0;
channel: sequence, l=3.0, refer=entry;
b, at=0.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
    lattice.set_all_string_attribute("extractor_type", "libff")
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


##########################################################################################


def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        ref_part, macroparticles, realparticles
    )
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    # populate 4x4 grid
    i = 0
    for x in np.linspace(-xoffs, xoffs, 5):
        for y in np.linspace(-yoffs, yoffs, 5):
            lp[i, 0] = x
            lp[i, 1] = xpoffs
            lp[i, 2] = y
            i = i + 1
    bunch.checkin_particles()
    return sim


##########################################################################################


def test_edge1(prop_fixture_edgekick):

    sim1 = create_simulator(
        prop_fixture_edgekick.get_lattice().get_reference_particle()
    )
    bunch1 = sim1.get_bunch()

    bunch1.checkout_particles()
    lp1 = bunch1.get_particles_numpy()
    lpinit = lp1.copy()
    # for ix in range(5):
    #     for iy in range(5):
    #         print(f'{ix}, {iy}: {lpinit[5*ix+iy, 0]}, {lpinit[5*ix+iy, 2]}')

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, True
    )
    prop_fixture_edgekick.propagate(sim1, simlog, nturns)

    bunch1.checkout_particles()
    lp1 = bunch1.get_particles_numpy()

    # particles started at y = 0 should have no py with edge kicks
    # all the particles with y != 0 should have nonzero py with edge kicks
    # index is ix*5+iy
    # y = 0 is iy = 2
    for ix in range(5):
        for iy in range(5):
            # print(f'{ix}, {iy} ({lpinit[5*ix+iy, 0]}, {lpinit[5*ix+iy, 2]}): lp1[{5*ix+iy}, 3]: ', lp1[5*ix+iy, 3])
            if iy == 2:
                assert lp1[5 * ix + 2, 3] == pytest.approx(0.0, abs=1.0e-15)
            else:
                assert abs(lp1[5 * ix + iy, 3]) > 1.0e-12


##########################################################################################


def test_edge2(prop_fixture_nokick):

    sim1 = create_simulator(prop_fixture_nokick.get_lattice().get_reference_particle())
    bunch1 = sim1.get_bunch()

    bunch1.checkout_particles()
    lp1 = bunch1.get_particles_numpy()
    lpinit = lp1.copy()
    # for ix in range(5):
    #     for iy in range(5):
    #         print(f'{ix}, {iy}: {lpinit[5*ix+iy, 0]}, {lpinit[5*ix+iy, 2]}')

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, True
    )
    prop_fixture_nokick.propagate(sim1, simlog, nturns)

    bunch1.checkout_particles()
    lp1 = bunch1.get_particles_numpy()

    # with no edge effects, there should be no particles with significant py
    for ix in range(5):
        for iy in range(5):
            # print(f'{ix}, {iy} ({lpinit[5*ix+iy, 0]}, {lpinit[5*ix+iy, 2]}): lp1[{5*ix+iy}, 3]: ', lp1[5*ix+iy, 3])
            assert abs(lp1[5 * ix + iy, 3]) == pytest.approx(0.0, abs=1.0e-15)


##########################################################################################
### it seems like the current implementation of bends does not perform edge kicks$%^&&*!!


def main():
    pf1 = prop_fixture_edgekick()
    pf2 = prop_fixture_nokick()
    test_edge1(pf1)
    test_edge2(pf2)


if __name__ == "__main__":
    main()
