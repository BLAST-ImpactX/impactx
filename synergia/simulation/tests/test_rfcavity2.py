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
    print("lattice: ")
    print(prop_fixture.get_lattice())

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

    # CDT should reflect the new velocity
    new_beta = new_p / new_E
    L = prop_fixture.get_lattice().get_length()
    assert lp[0, 4] == pytest.approx(L * (1 / new_beta - 1 / orig_beta))
    print("lp[0]: ", lp[0, :])

    print(
        "bunch design energy: ",
        bunch.get_design_reference_particle().get_total_energy(),
    )
    print("bunch energy: ", bunch.get_reference_particle().get_total_energy())

    # check new momenta
    mom1_i = (1 + dpop_offset) * orig_p
    e1_i = np.sqrt(mom1_i**2 + mp**2)
    e1_f = e1_i + expected_delta_E
    mom1_f_should_be = np.sqrt(e1_f**2 - mp**2)
    # check it
    assert mom1_f_should_be == pytest.approx((1 + lp[1, 5]) * new_p)

    mom2_i = (1 - dpop_offset) * orig_p
    e2_i = np.sqrt(mom2_i**2 + mp**2)
    e2_f = e2_i + expected_delta_E
    mom2_f_should_be = np.sqrt(e2_f**2 - mp**2)
    # check it
    assert mom2_f_should_be == pytest.approx((1 + lp[2, 5]) * new_p)


def main():
    pf = prop_fixture()
    test_accel1(pf)


if __name__ == "__main__":
    main()
