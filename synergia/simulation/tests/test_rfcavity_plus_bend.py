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
b: sbend, angle=pi/2, l=pi/2;
channel: sequence, l=pi/2, refer=entry;
rfc, at=0.0;
b, at=0.0;
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
    assert lp[0, 5] == 0.0

    # check offset of increased momentum particle
    # Need the old and new radius of curvature for the sbend
    # With new momentum, radius changes by D=newR-oldR, the x offset after
    # 90 degree bend is R*sqrt(1 + 2*D/R) - R
    # can also be expressed as R*sqrt(2*newR/oldR-1)-R

    lattice = prop_fixture.get_lattice()
    bend = None
    for elem in lattice.get_elements():
        if elem.get_type() == synergia.lattice.element_type.sbend:
            bend = elem
            break
    if not bend:
        raise RuntimeError("no bend element found")

    # L = R theta => R = L/theta
    oldR = bend.get_length() / bend.get_bend_angle()

    # newR/oldR = new-p/old-p

    newR = oldR * new_p / orig_p
    print("oldR: ", oldR, " newR: ", newR, " difference: ", newR - oldR)

    x = oldR * (np.sqrt(2 * (newR / oldR) - 1)) - oldR
    print("x offset: ", x)
    assert lp[0, 0] == pytest.approx(x)

    # check cdt
    new_bend_angle = np.arccos((newR - oldR) / newR)
    print("new bend angle: ", new_bend_angle)
    new_beta = bunch.get_reference_particle().get_beta()
    old_beta = bunch.get_design_reference_particle().get_beta()
    old_ctime = bend.get_length() / old_beta
    new_ctime = newR * new_bend_angle / new_beta
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
