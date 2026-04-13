#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

if hasattr(synergia.bunch.Bunch, "read_openpmd_file"):
    import openpmd_api as io
else:
    import h5py

macroparticles = 16
realparticles = 4.0e10
# lag 1/120 is a phase angle of 2pi/120 or pi/60 or 3 degrees
# V = 0.2 MV * sin(pi/60) =
turn_voltage = 1.0e-3  # 1.0e-3 GV/turn
expected_delta_E = turn_voltage * np.sin(np.pi / 60)
print("expected delta E/turn: ", expected_delta_E)
nturns = 10


# prop_fixture is a propagator
# can't run this test with the fixture because I need the propagator
# and simulator to be deleted so I can read the contents of the
# generated hdf5 files
#
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


def prop_accel_bunch(prop_fixture):
    refpart = prop_fixture.get_lattice().get_reference_particle()
    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())
    diag_part = synergia.bunch.Diagnostics_particles("diag_part.h5")
    sim.reg_diag_per_turn(diag_part)

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

    return context.energies


def check_energies(energies):
    if hasattr(synergia.bunch.Bunch, "read_openpmd_file"):
        # test if openpmd

        s = io.Series("diag_part.h5", io.Access_Type.read_only)
        iters = s.iterations
        N = len(iters)
        print(N, " iterations in file")
        for k in range(N):
            mass = iters[k].particles["bunch_particles"].get_attribute("mass")
            # file mass is stored in kg
            # strangely enough, if I use the kg_to_GeV constant in this calculation instead
            # of explicitly converting using c**2/e then the energy differs in the last
            # decimal place.
            mass = (
                mass
                * 1.0e-9
                * synergia.foundation.pconstants.c**2
                / synergia.foundation.pconstants.e
            )
            # mass = mass * synergia.foundation.pconstants.kg_to_GeV
            gamma_ref = iters[k].particles["bunch_particles"].get_attribute("gamma_ref")
            step_energy = mass * gamma_ref
            print(
                "iteration",
                k,
                " gamma_ref: ",
                gamma_ref,
                ", step energy: ",
                step_energy,
                ", prop energy: ",
                energies[k],
                flush=True,
            )
            assert energies[k] == pytest.approx(step_energy, rel=1.0e-12)
        s.close()
        os.remove("diag_part.h5")

    else:
        # legacy branch
        N = len(energies)
        for k in range(N):
            h5 = h5py.File(f"diag_part_{k:05d}.h5", "r")
            mass = h5.get("mass")[()]
            step_pz = h5.get("pz")[()]
            step_energy = np.sqrt(step_pz**2 + mass**2)
            print(
                "interation",
                k,
                "step_energy: ",
                step_energy,
                ", prop_energy: ",
                energies[k],
            )
            assert step_energy == pytest.approx(energies[k], rel=1.0e-12)
            h5.close()
            os.remove(f"diag_part_{k:05d}.h5")

    return


def test_accel_write_openpmd():
    pf = prop_fixture()
    energies = prop_accel_bunch(pf)
    del pf
    check_energies(energies)


def main():
    test_accel_write_openpmd()


if __name__ == "__main__":
    main()
