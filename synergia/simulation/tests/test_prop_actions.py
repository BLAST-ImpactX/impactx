#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10
# lag 1/120 is a phase angle of 2pi/120 or pi/60 or 3 degrees
# V = 0.2 MV * sin(pi/60) =
expected_delta_E = 0.0002 * np.sin(np.pi / 60)
print("expected delta E/turn: ", expected_delta_E)
nturns = 1


@pytest.fixture
def prop_fixture():
    fodo_madx = """
beam, particle=proton,pc=0.75*pmass;

f: quadrupole, l=1.0, k1=0.0625;
d: quadrupole, l=1.0, k1=-0.0625;
rfc: rfcavity, l=0.0, volt=0.2, harmon=1, lag=(1/120.0);
m1: marker;
m2: marker;

fodo: sequence, l=20.0, refer=centre;
fodo_0: m1, at=0.0;
fodo_1: f, at=1.0;
fodo_2: d, at=9.0;
m2, at=10.0;
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
    lattice.set_all_string_attribute("extractor_type", "libff")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(
        ref_part, macroparticles, realparticles
    )
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch.checkout_particles()
    return sim


def create_propagator(lattice):
    sc_ops = synergia.collective.Dummy_CO_options()

    stepper = synergia.simulation.Split_operator_stepper_elements(sc_ops, 1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


def test_modify_lattice(prop_fixture):
    propagator = prop_fixture
    lattice = propagator.get_lattice()

    sim = create_simulator(lattice.get_reference_particle())

    # turn and action method
    def turn_end_action(sim, lattice, turn):
        lattice.get_elements()[0].set_double_attribute("foo", 3.125 * (turn + 1))
        print("test_modify lattice turn ", turn)

    # end of turn end action method

    sim.reg_prop_action_turn_end(turn_end_action)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    propagator.propagate(sim, simlog, 10)
    print("foo attribute: ", lattice.get_elements()[0].get_double_attribute("foo"))
    print(
        "foo attribute through propagator: ",
        propagator.get_lattice().get_elements()[0].get_double_attribute("foo"),
    )
    assert (
        propagator.get_lattice().get_elements()[0].get_double_attribute("foo") == 31.25
    )
    # assert False


def test_modify_lattice_energy(prop_fixture):
    propagator = prop_fixture
    lattice = propagator.get_lattice()

    sim = create_simulator(lattice.get_reference_particle())

    # turn and action method
    def turn_end_action(sim, lattice, turn):
        lattice.set_lattice_energy(8.0 + (turn + 1) * 0.25)
        print("test_modify lattice_energy turn ", turn)

    # end of turn end action method

    sim.reg_prop_action_turn_end(turn_end_action)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    propagator.propagate(sim, simlog, 1)
    print("lattice energy: ", propagator.get_lattice().get_lattice_energy())

    assert propagator.get_lattice().get_lattice_energy() == 8.25
    # assert False


# this test activates a context class that can maintain separate state but is accessible within the
# prop_action. This test will increment a counter for odd numbered turns.
def test_context(prop_fixture):
    propagator = prop_fixture
    lattice = propagator.get_lattice()

    sim = create_simulator(lattice.get_reference_particle())

    class context:
        odd_turn_count = 0

    # turn and action method
    def turn_end_action(sim, lattice, turn):
        if turn % 2 == 1:
            context.odd_turn_count = context.odd_turn_count + 1
        print(
            "test_context: turn  ", turn, ", odd_turn_count: ", context.odd_turn_count
        )

    # end of turn end action method

    sim.reg_prop_action_turn_end(turn_end_action)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    propagator.propagate(sim, simlog, 5)

    print("after propagate: odd_turn_count: ", context.odd_turn_count)
    assert context.odd_turn_count == 2
    # assert False


# --------------------------------------------------------------------------------------------------------
def fixmarker(inlattice, turn):
    print("lattice id: ", id(inlattice))
    for elem in inlattice.get_elements():
        if elem.get_name() == "m2":
            elem.set_double_attribute("foo", 2.5 * (turn + 1))


def test_modify_lattice2(prop_fixture):
    propagator = prop_fixture
    lattice = propagator.get_lattice()

    sim = create_simulator(lattice.get_reference_particle())

    # turn and action method
    def turn_end_action(sim, lattice, turn):
        print("turn_end_action, turn: ", turn, ", lattice: ", id(lattice))
        fixmarker(lattice, turn)

    # end of turn end action method

    sim.reg_prop_action_turn_end(turn_end_action)

    print("propagator.lattice id: ", id(propagator.get_lattice()))
    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    propagator.propagate(sim, simlog, 10)
    for elem in propagator.get_lattice().get_elements():
        if elem.get_name() == "m2":
            assert elem.get_double_attribute("foo") == 25.0
    # assert False


# --------------------------------------------------------------------------------------------------------
def main():
    pf = prop_fixture()
    test_modify_lattice2(pf)


if __name__ == "__main__":
    main()
