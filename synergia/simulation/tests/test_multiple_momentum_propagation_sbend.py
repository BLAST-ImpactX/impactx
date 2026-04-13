#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles = 16
realparticles = 4.0e10

nturns = 1

mp = synergia.foundation.pconstants.mp


# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    channel_madx = """
beam, particle=proton,energy=8.0+pmass;
b: sbend, angle=pi/2, l=pi/2;
channel: sequence, l=pi/2, refer=entry;
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


# two trains start out with the same reference particle and design reference momentum
# because they are propagating in the
# the same accelerator, but the secondary bunch will have it bunch energy increase
# by 100 MeV to mimic slip stacking. The bunches propagate in the same pipe but one
# secondary moves faster than the primary so its cdt will go negative compared to
# primary.
def create_simulator(ref_part_pri, ref_part_sec):
    sim = synergia.simulation.Bunch_simulator.create_two_trains_simulator(
        ref_part_pri, ref_part_sec, macroparticles, realparticles
    )
    bunch_pri = sim.get_bunch(0, 0)
    bunch_pri.checkout_particles()
    lp = bunch_pri.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch_pri.checkin_particles()

    design_energy = ref_part_pri.get_total_energy()
    bunch_sec = sim.get_bunch(1, 0)
    bunch_sec.get_reference_particle().set_total_energy(design_energy + 0.100)
    bunch_sec.checkout_particles()
    lp = bunch_sec.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch_sec.checkin_particles()

    return sim


def test_co_propagate(prop_fixture):
    # print('lattice: ')
    # print(prop_fixture.get_lattice())

    refpart = prop_fixture.get_lattice().get_reference_particle()

    sim = create_simulator(refpart, refpart)

    simlog = synergia.utils.parallel_utils.Logger(
        0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False
    )
    prop_fixture.propagate(sim, simlog, nturns)

    bunch_pri = sim.get_bunch(0, 0)
    bunch_sec = sim.get_bunch(1, 0)

    # primary bunch should be locked to the main energy so its cdt will be 0
    # the secondary bunch will slip compared to primary so its particle's cdt
    # will be negative.

    lattice = prop_fixture.get_lattice()
    L = lattice.get_length()

    beta_pri = bunch_pri.get_reference_particle().get_beta()
    beta_sec = bunch_sec.get_reference_particle().get_beta()
    pref_pri = bunch_pri.get_reference_particle().get_momentum()
    pref_sec = bunch_sec.get_reference_particle().get_momentum()
    design_pref_pri = bunch_pri.get_design_reference_particle().get_momentum()
    design_pref_sec = bunch_sec.get_design_reference_particle().get_momentum()

    print("beta_pri: ", beta_pri)
    print("beta_sec: ", beta_sec)
    print("pref_pri: ", pref_pri)
    print("design pref pri: ", design_pref_pri)
    print("pref_sec: ", pref_sec)
    print("design_pref_sec: ", design_pref_sec)

    bunch_pri.checkout_particles()
    lp_pri = bunch_pri.get_particles_numpy()
    bunch_sec.checkout_particles()
    lp_sec = bunch_sec.get_particles_numpy()

    assert lp_pri[0, 4] == pytest.approx(0.0)

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

    newR = oldR * pref_sec / pref_pri
    print("oldR: ", oldR, " newR: ", newR, " difference: ", newR - oldR)

    x = oldR * (np.sqrt(2 * (newR / oldR) - 1)) - oldR
    assert lp_sec[0, 0] == pytest.approx(x)

    new_bend_angle = np.arccos((newR - oldR) / newR)
    new_length = newR * new_bend_angle

    ctime_diff = new_length / beta_sec - bend.get_length() / beta_pri
    print("after propagation cdt: ", lp_sec[0, 4])
    print("calculated cdt: ", ctime_diff)

    assert ctime_diff == pytest.approx(lp_sec[0, 4])

    # assert False


def main():
    pf = prop_fixture()
    test_accel1(pf)


if __name__ == "__main__":
    main()
