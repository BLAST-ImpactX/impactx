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
rfc: rfcavity, l=0.0, volt=0.2, harmon=1, lag=0;
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
    # calculate dp/p for additional energy 0.0001
    origE = bunch.get_design_reference_particle().get_total_energy()
    origp = bunch.get_design_reference_particle().get_momentum()
    newE = origE + expected_delta_E
    newp = np.sqrt(newE**2 - mp**2)
    dpop = newp / origp - 1
    lp[1, 5] = dpop
    bunch.checkin_particles()
    return sim


def test_accel1(prop_fixture):
    # print('lattice: ')
    # print(prop_fixture.get_lattice())

    refpart = prop_fixture.get_lattice().get_reference_particle()

    origE = refpart.get_total_energy()
    origp = refpart.get_momentum()
    newE = origE + expected_delta_E
    newp = np.sqrt(newE**2 - mp**2)
    dpop = newp / origp - 1

    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())
    bunch = sim.get_bunch()

    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

    print("dpop of test particle: ", lp[1, 5])

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    prop_fixture.propagate(sim, simlog, nturns)

    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

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

    newR = oldR * newp / origp
    print("oldR: ", oldR, " newR: ", newR, " difference: ", newR - oldR)

    x = oldR * (np.sqrt(2 * (newR / oldR) - 1)) - oldR
    assert lp[1, 0] == pytest.approx(x)

    new_bend_angle = np.arccos((newR - oldR) / newR)
    new_length = newR * new_bend_angle
    oldbeta = refpart.get_beta()
    newbeta = newp / newE
    ctime_diff = new_length / newbeta - bend.get_length() / oldbeta
    print("after propagation cdt: ", lp[1, 4])
    print("calculated cdt: ", ctime_diff)

    assert ctime_diff == pytest.approx(lp[1, 4])


def main():
    pf = prop_fixture()
    test_accel1(pf)


if __name__ == "__main__":
    main()
