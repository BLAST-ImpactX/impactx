#!/usr/bin/env python

import numpy as np
import pytest

import synergia

mass = 1.0
betagamma = 3 / 4
gamma = 5 / 4
energy = mass * gamma

macroparticles = 16
realparticles = 5e10


@pytest.fixture
def bunch_fixture():
    refpart = synergia.foundation.Reference_particle(1, mass, energy)
    bunch = synergia.bunch.Bunch(
        refpart, macroparticles, realparticles, synergia.utils.Commxx()
    )
    bunch.checkin_particles()
    return bunch


def test_mask_access(bunch_fixture):
    bunch_fixture.checkout_particles()
    masks = bunch_fixture.get_particle_masks_numpy()


def test_mask_initial_state(bunch_fixture):
    # all masks should be 1 when bunch is created
    bunch_fixture.checkout_particles()
    masks = bunch_fixture.get_particle_masks_numpy()
    assert masks.sum() == macroparticles


def test_mask_after_cuts():
    refpart = synergia.foundation.Reference_particle(1, mass, energy)
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        refpart, macroparticles, realparticles, synergia.utils.Commxx()
    )
    bunch = sim.get_bunch()
    lp = bunch.get_particles_numpy()
    # populate between -2mm to +2mm

    lp[:, 0] = 0.002 * (np.arange(macroparticles, dtype="d") / macroparticles - 0.5)
    bunch.checkin_particles()

    # create a lattice with an aperture
    lattice = synergia.lattice.Lattice("foo")
    m = synergia.lattice.Lattice_element("marker", "m")
    m.set_string_attribute("aperture_type", "rectangular")
    m.set_double_attribute("rectangular_aperture_width", 0.00030 * 2)
    m.set_double_attribute("rectangular_aperture_height", 0.002)
    lattice.append(m)
    lattice.set_reference_particle(refpart)

    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, True
    )

    propagator.propagate(sim, simlog, 1)
    bunch.checkout_particles()
    masks = bunch.get_particle_masks_numpy()
    assert masks.sum() == 5


# test access mask
def run_test1():
    bf = bunch_fixture()
    test_mask_access(bf)


# test initial state of mask
def run_test2():
    bf = bunch_fixture()
    test_mask_initial_state(bf)


# test masks after cutting some particles in aperture
def run_test3():
    test_mask_after_cuts()


# for running from the command line, comment out pytest.fixture
if __name__ == "__main__":
    run_test1()
    run_test2()
    run_test3()
