#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest
import h5py

macroparticles = 16
realparticles = 4.0e10
# lag 1/120 is a phase angle of 2pi/120 or pi/60 or 3 degrees
# V = 0.2 MV * sin(pi/60) =
expected_delta_E = 0.0002 * np.sin(np.pi / 60)
print("expected delta E/turn: ", expected_delta_E)
nturns = 10
mp = synergia.foundation.pconstants.mp


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    fodo_madx = """
beam, particle=proton,pc=0.75*pmass;

f: quadrupole, l=1.0, k1=0.0625;
d: quadrupole, l=1.0, k1=-0.0625;
rfc: rfcavity, l=0.0, volt=0.2, harmon=1, lag=(1/120.0);

fodo: sequence, l=20.0, refer=centre;
fodo_1: f, at=1.0;
fodo_2: d, at=9.0;
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


def test_lattice_energy(prop_fixture):
    energy = prop_fixture.get_lattice().get_lattice_energy()
    assert energy == pytest.approx(synergia.foundation.pconstants.mp * 1.25)


def test_lattice_length(prop_fixture):
    assert prop_fixture.get_lattice().get_length() == 20.0


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


def test_accel1(prop_fixture):
    # abort test if openPMD IO is active
    names = dir(synergia.bunch.Bunch)
    if names.count("write_openpmd_file"):
        return True

    refpart = prop_fixture.get_lattice().get_reference_particle()
    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())

    # E^2 = p^2 + m^2
    # 2 E dE = 2 p dp
    # expected dp = E/p dE
    expected_dp = refpart.get_total_energy() / refpart.get_momentum() * expected_delta_E

    lattice = prop_fixture.get_lattice()

    Elat0 = lattice.get_lattice_energy()
    Ebun0 = sim.get_bunch().get_design_reference_particle().get_total_energy()
    assert Elat0 == pytest.approx(Ebun0, 1.0e-10)

    # diagnostics
    diag_full2 = synergia.bunch.Diagnostics_full2("diag.h5")
    sim.reg_diag_per_turn(diag_full2)

    diag_bt = synergia.bunch.Diagnostics_bulk_track("tracks.h5", 4, 0)
    sim.reg_diag_per_turn(diag_bt)

    diag_part = synergia.bunch.Diagnostics_particles("particles.h5", 4)
    sim.reg_diag_per_turn(diag_part)

    # turn and action method
    def turn_end_action(sim, lattice, turn):
        bunch = sim.get_bunch()
        bunch_design_E = bunch.get_design_reference_particle().get_total_energy()
        bunch_E = bunch.get_reference_particle().get_total_energy()
        lattice_E = lattice.get_lattice_energy()

        print("turn_end_action: enter: bunch_design_E: ", bunch_design_E)
        print("turn_end_action: enter: lattice_E: ", lattice_E)
        print("turn_end_action: enter: bunch_E: ", bunch_E)

        # after RF cavity, the bunch energy should have increased but
        # neither the bunch design energy nor the lattice energy
        # will have increased.

        # set the bunch design energy and lattice energy to match the bunch
        # energy
        bunch.get_design_reference_particle().set_total_energy(bunch_E)
        lattice.set_lattice_energy(bunch_E)

        print(
            "turn_end_action: exit: bunch_design_E: ",
            bunch.get_design_reference_particle().get_total_energy(),
        )
        print(
            "turn_end_action: exit: lattice_E: ",
            lattice.get_reference_particle().get_total_energy(),
        )
        print(
            "turn_end_action: exit: bunch_E: ",
            bunch.get_reference_particle().get_total_energy(),
        )

        # tune lattice
        synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)

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

    # check pz in the diagnostics output
    del sim
    del prop_fixture

    # calculate the momentas
    energies = np.zeros(nturns + 1)
    energies[0] = Ebun0
    for i in range(1, nturns + 1):
        energies[i] = energies[i - 1] + expected_delta_E
    pzs = np.sqrt(energies**2 - mp**2)

    h5diag = h5py.File("diag.h5", "r")
    # are there the correct number of entries?
    print("number of pz entries in diag: ", h5diag.get("pz").shape[0])
    assert h5diag.get("pz").shape[0] == nturns + 1
    # Is the step between entries correct?
    print("pz[1]: ", h5diag.get("pz")[1])
    print("pz[0]: ", h5diag.get("pz")[0])
    print("expected dp: ", expected_dp)
    # This test only works as-is with the first step because expected_dp changes
    assert h5diag.get("pz")[1] - h5diag.get("pz")[0] == pytest.approx(
        expected_dp, rel=1.0e-5
    )

    for i in range(nturns + 1):
        assert h5diag.get("pz")[i] == pytest.approx(pzs[i])

    # check pz in tracks file
    h5trks = h5py.File("tracks.h5", "r")
    assert h5trks.get("track_pz").shape[0] == nturns + 1
    assert h5trks.get("track_pz")[1] - h5trks.get("track_pz")[0] == pytest.approx(
        expected_dp, rel=1.0e-5
    )
    for i in range(nturns + 1):
        assert h5trks.get("track_pz")[i] == pytest.approx(h5diag.get("pz")[i])

    # test pz in particles file
    for i in range(nturns + 1):
        h5parts = h5py.File(f"particles_{i:05d}.h5", "r")
        assert h5trks.get("track_pz")[i] == h5parts.get("pz")
        h5parts.close()

    h5diag.close()
    h5trks.close()


def main():
    pf = prop_fixture()
    test_accel1(pf)


if __name__ == "__main__":
    main()
