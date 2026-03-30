#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10
# lag 1/120 is a phase angle of 2pi/120 or pi/60 or 3 degrees
# V = 0.2 MV * sin(pi/60) =
turn_voltage = 1.0e-3  # 1.0e-3 GV/turn
expected_delta_E = turn_voltage * np.sin(np.pi / 60)
print("expected delta E/turn: ", expected_delta_E)
nturns = 100


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    # booster-like lattice
    booster_madx = """
ncells=24;
turn_voltage=1.0; ! 1 MV /turn
beam, particle=proton,energy=pmass+0.8;

f: sbend, l=2.0, angle=(pi/(2*ncells)), k1=1/16.2;
d: sbend, l=2.0, angle=(pi/(2*ncells)), k1=-1/16.7;
!f: quadrupole, l=2.0, k1=0.0625;
!d: quadrupole, l=2.0, k1=-0.0625;
rfc: rfcavity, l=0.0, volt=turn_voltage/ncells, harmon=96, lag=(1/120.0);


cell: sequence, l=20.0, refer=centre;
fodo_1: f, at=1.5;
fodo_2: d, at=8.5;
fodo_3: d, at=11.5;
fodo_4: f, at=18.5;
fodo_5: rfc, at=20.0;
endsequence;

booster: sequence, l=480.0, refer=entry;
cell, at=0.0;
cell, at=20.0;
cell, at=40.0;
cell, at=60.0;
cell, at=80.0;
cell, at=100.0;
cell, at=120.0;
cell, at=140.0;
cell, at=160.0;
cell, at=180.0;
cell, at=200.0;
cell, at=220.0;
cell, at=240.0;
cell, at=260.0;
cell, at=280.0;
cell, at=300.0;
cell, at=320.0;
cell, at=340.0;
cell, at=360.0;
cell, at=380.0;
cell, at=400.0;
cell, at=420.0;
cell, at=440.0;
cell, at=460.0;
endsequence;

"""

    reader = synergia.lattice.MadX_reader()
    reader.parse(booster_madx)
    lattice = reader.get_lattice("booster")
    lattice.set_all_string_attribute("extractor_type", "libff")
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator


def test_lattice_energy(prop_fixture):
    energy = prop_fixture.get_lattice().get_lattice_energy()
    assert energy == pytest.approx(synergia.foundation.pconstants.mp + 0.8)


def test_lattice_length(prop_fixture):
    assert prop_fixture.get_lattice().get_length() == 480.0


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


def test_accel2(prop_fixture):
    refpart = prop_fixture.get_lattice().get_reference_particle()
    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())

    lattice = prop_fixture.get_lattice()

    Elat0 = lattice.get_lattice_energy()
    Ebun0 = sim.get_bunch().get_design_reference_particle().get_total_energy()
    assert Elat0 == pytest.approx(Ebun0, 1.0e-10)

    class context:
        max_dpop = 0.0
        max_cdt = 0.0
        energies = [Elat0]  # save lattice energy at each turn

    # turn and action method
    def turn_end_action(sim, lattice_in, turn):
        bunch = sim.get_bunch()
        bunch_design_E = bunch.get_design_reference_particle().get_total_energy()
        bunch_E = bunch.get_reference_particle().get_total_energy()
        lattice_E = lattice_in.get_lattice_energy()

        # print('turn_end_action: enter: bunch_design_E: ', bunch_design_E)
        # print('turn_end_action: enter: lattice_E: ', lattice_E)
        # print('turn_end_action: enter: bunch_E: ', bunch_E)

        # after RF cavity, the bunch energy should have increased but
        # neither the bunch design energy nor the lattice energy
        # will have increased.

        # set the bunch design energy and lattice energy to match the bunch
        # energy
        bunch.get_design_reference_particle().set_total_energy(bunch_E)
        lattice_in.set_lattice_energy(bunch_E)

        # save new lattice energy
        context.energies.append(lattice_in.get_lattice_energy())
        # print('turn_end_action: exit: bunch_design_E: ', bunch.get_design_reference_particle().get_total_energy())
        # print('turn_end_action: exit: lattice_E: ', lattice.get_reference_particle().get_total_energy())
        # print('turn_end_action: exit: bunch_E: ', bunch.get_reference_particle().get_total_energy())

        # tune lattice
        synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice_in)

        # check frequency matches new energy
        beta1 = lattice_in.get_reference_particle().get_beta()
        beta2 = bunch.get_reference_particle().get_beta()
        assert beta1 == pytest.approx(beta2)
        freq = 96 * beta1 * synergia.foundation.pconstants.c / lattice_in.get_length()
        assert freq == pytest.approx(
            synergia.simulation.Lattice_simulator.get_rf_frequency(lattice)
        )

        # The central particles should stay close to 0 in energy and time
        bunch.checkout_particles()
        lp = bunch.get_particles_numpy()

        if abs(lp[0, 4]) > context.max_cdt:
            context.max_cdt = lp[0, 4]
        if abs(lp[0, 5]) > context.max_dpop:
            context.max_dpop = lp[0, 5]

    # end of turn end action method

    sim.reg_prop_action_turn_end(turn_end_action)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    prop_fixture.propagate(sim, simlog, nturns)

    Ebun1 = sim.get_bunch().get_reference_particle().get_total_energy()
    Elat1 = prop_fixture.get_lattice().get_lattice_energy()
    print("(Ebun1-Ebun0)/expected_delta_E: ", (Ebun1 - Ebun0) / expected_delta_E)
    print("(Elat1-Elat0)/expected_delta_E: ", (Elat1 - Elat0) / expected_delta_E)
    assert (Ebun1 - Ebun0) / expected_delta_E == pytest.approx(nturns)
    assert (Elat1 - Elat0) / expected_delta_E == pytest.approx(nturns)

    assert context.max_cdt < 1.0e-2
    assert context.max_dpop < 1.0e-5

    # check saved energy gain each turn
    for i in range(1, len(context.energies)):
        assert (context.energies[i] - context.energies[i - 1]) == pytest.approx(
            expected_delta_E
        )


def main():
    pf = prop_fixture()
    test_accel2(pf)


if __name__ == "__main__":
    main()
