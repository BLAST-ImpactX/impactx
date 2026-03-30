#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10
# lag 1/12 is a phase angle of 2pi/12 or pi/6 or 30 degrees
# V = 0.2 MV * sin(pi/6) =
expected_delta_E = 0.0002 * np.sin(np.pi / 6)
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
rfc: rfcavity, l=0.0, volt=0.2, harmon=1, lag=(1/12.0);
channel: sequence, l=20.0, refer=centre;
rfc, at=0.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice("channel")
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
    lp[1, 5] = dpop_offset  # particle with dp/p offset
    lp[2, 5] = -dpop_offset  # particle with negative dp/p offset
    lp[3, 1] = transmom_offset  # particle with transverse x momentum
    lp[4, 3] = -transmom_offset  # particle with transverse y momentum

    lp[5, 1] = dpop_offset  # particle with both trans xmomentum and dp/p
    lp[5, 5] = transmom_offset

    lp[6, 3] = -transmom_offset  # particle with both trans y momentum and negative dp/p
    lp[6, 5] = -dpop_offset

    bunch.checkin_particles()
    return sim


def test_accel1(prop_fixture):
    # print('lattice: ')
    # print(prop_fixture.get_lattice())

    refpart = prop_fixture.get_lattice().get_reference_particle()

    orig_E = refpart.get_total_energy()
    orig_p = refpart.get_momentum()
    orig_beta = orig_p / orig_E
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

    # after acceleration, the dp/p of the particles should still be 0
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    assert lp[0, 5] == pytest.approx(0.0)

    # CDT should reflect the new velocity
    new_E = sim.get_bunch().get_reference_particle().get_total_energy()
    new_p = sim.get_bunch().get_reference_particle().get_momentum()
    print("new_E: ", new_E)
    print("new_p: ", new_p)
    assert new_E - orig_E == pytest.approx(nturns * expected_delta_E)

    new_beta = new_p / new_E
    L = prop_fixture.get_lattice().get_length()
    assert lp[0, 4] == pytest.approx(L * (1 / new_beta - 1 / orig_beta))
    print("lp[0]: ", lp[0, :])

    print(
        "bunch design energy: ",
        bunch.get_design_reference_particle().get_total_energy(),
    )
    print("bunch energy: ", bunch.get_reference_particle().get_total_energy())

    # check new momentum for particle 1 with dpop offset
    mom1_i = (1 + dpop_offset) * orig_p  # original momentum
    e1_i = np.sqrt(mom1_i**2 + mp**2)  # original energy
    e1_f = e1_i + expected_delta_E  # new energy after accel
    mom1_f_should_be = np.sqrt(e1_f**2 - mp**2)  # new momentum after accel
    # WTF why doesn't e1_f printout show energy gain?
    assert expected_delta_E != pytest.approx(0.0)
    assert e1_f == pytest.approx(e1_i + expected_delta_E)

    # printout for debugging test on v100 GPU
    print("dpop_offset: ", dpop_offset)
    print("mom1_i: ", mom1_i)
    print("e1_i: ", e1_i)
    print("e1_f: ", e1_f)
    print("mom1_f_should_be: ", mom1_f_should_be)
    print("lp[1, 5]: ", lp[1, 5])
    print("mom1_f_should_be approx: ", (1 + lp[1, 5]) * new_p)

    # check it
    assert mom1_f_should_be == pytest.approx((1 + lp[1, 5]) * new_p)

    # check new momentum for particle with negative dpop offset
    mom2_i = (1 - dpop_offset) * orig_p  # original momentum
    e2_i = np.sqrt(mom2_i**2 + mp**2)  # original energy
    e2_f = e2_i + expected_delta_E  # energy after accel
    mom2_f_should_be = np.sqrt(e2_f**2 - mp**2)  # momentum after accel
    # check it
    assert mom2_f_should_be == pytest.approx((1 + lp[2, 5]) * new_p)

    # check transverse momenta properly transformed
    assert lp[3, 1] == pytest.approx(transmom_offset * orig_p / new_p)
    # total momentum is design momentum so it should continue to be design momentum
    assert new_p == pytest.approx((1 + lp[3, 5]) * new_p)

    assert lp[4, 3] == pytest.approx(-transmom_offset * orig_p / new_p)
    # total momentum is design momentum so it should continue to be design momentum
    assert new_p == pytest.approx((1 + lp[4, 5]) * new_p)

    # check particles with both trans and long momentum
    assert lp[5, 1] == pytest.approx(transmom_offset * orig_p / new_p)
    assert mom1_f_should_be == pytest.approx((1 + lp[5, 5]) * new_p)

    assert lp[6, 3] == pytest.approx(-transmom_offset * orig_p / new_p)
    assert mom2_f_should_be == pytest.approx((1 + lp[6, 5]) * new_p)


def main():
    pf = prop_fixture()
    test_accel1(pf)


if __name__ == "__main__":
    main()
