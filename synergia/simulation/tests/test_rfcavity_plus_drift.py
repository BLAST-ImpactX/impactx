#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10
# lag 1/12 is a phase angle of 2pi/120 or pi/60 or 3 degrees
# V = 0.2 MV * sin(pi/60) =
expected_delta_E = 0.0002 * np.sin(np.pi / 60)
print("expected delta E/turn: ", expected_delta_E)
nturns = 1
dpop_offset = 1.0e-3
transmom_offset = 0.001
mp = synergia.foundation.pconstants.mp


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    channel_madx = """
beam, particle=proton,pc=0.75*pmass;
rfc: rfcavity, l=0.0, volt=0.2, harmon=1, lag=(1/120.0);
channel: sequence, l=20.0, refer=entry;
rfc, at=0.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
    lattice.set_all_string_attribute("extractor_type", "libff")
    synergia.simulation.Lattice_simulator.tune_linear_lattice(lattice)
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


def test_accel1(prop_fixture):
    # print('lattice: ')
    # print(prop_fixture.get_lattice())

    refpart = prop_fixture.get_lattice().get_reference_particle()

    orig_E = refpart.get_total_energy()
    orig_p = refpart.get_momentum()
    print("orig_E: ", orig_E)
    print("orig_p: ", orig_p)
    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    bunch.checkin_particles()

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    prop_fixture.propagate(sim, simlog, nturns)

    new_E = sim.get_bunch().get_reference_particle().get_total_energy()
    new_p = sim.get_bunch().get_reference_particle().get_momentum()
    print("new_E: ", new_E)
    print("new_p: ", new_p)
    assert new_E - orig_E == pytest.approx(nturns * expected_delta_E)

    # after acceleration, the dp/p of the particles should still be 0
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    assert lp[0, 5] == pytest.approx(0.0)

    length = prop_fixture.get_lattice().get_length()

    # check cdt
    new_beta = bunch.get_reference_particle().get_beta()
    old_beta = bunch.get_design_reference_particle().get_beta()
    old_ctime = length / old_beta
    new_ctime = length / new_beta
    ctime_diff = new_ctime - old_ctime
    print("ctime diff: ", ctime_diff)
    print("lp[0, 4]: ", lp[0, 4])

    # The reference time for a particle should be based on
    # the design momentum of the bunch, not momentum after
    #  acceleration.
    assert ctime_diff == pytest.approx(lp[0, 4])

    # assert False


def main():
    pf = prop_fixture()
    test_accel1(pf)


if __name__ == "__main__":
    main()
