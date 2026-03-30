#!/usr/bin/env python

import sys
import os
import synergia
import synergia.tools
import numpy as np
import h5py
import pytest

#  just a little tester for the class
def test_three_bump():
    lattice = synergia.lattice.MadX_reader().get_lattice("model", "lattices/foborodobo128.madx")
    print("read lattice: ", len(lattice.get_elements()), " elements, length = ", lattice.get_length())
    hcorr_names = ('hc1', 'hc2', 'hc3')
    vcorr_names = ('vc1', 'vc2', 'vc3')
    three_bump = synergia.tools.Three_bump(lattice, 'm1', 'm2', hcorr_names, vcorr_names, 'm3', False)
    three_bump.information()

    elem_m3 = None
    for e in lattice.get_elements():
        if e.get_name() == 'm3':
            elem_m3 = e
            break
    if elem_m3 is None:
        raise RuntimeError('no m3 element found')
        
    target_x = 0.001
    target_y = -0.0005
    bump_settings = three_bump.set_bump((target_x, target_y))
    print("bump_settings: ", bump_settings[0], bump_settings[1], bump_settings[2], bump_settings[3], bump_settings[4], bump_settings[5])

    # propagate the whole lattice now
    comm = synergia.utils.Commxx()
    refpart = lattice.get_reference_particle()
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(refpart, 8, 0.5e11)
    bunch = sim.get_bunch(0, 0)
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch.checkin_particles()

    # register diagnostics
    diag_bump = synergia.bunch.Diagnostics_bulk_track("bump.h5", 1)
    diag_orbit = synergia.bunch.Diagnostics_bulk_track("fulllattice.h5", 1)
    sim.reg_diag_at_element(diag_bump, elem_m3)
    sim.reg_diag_per_turn(diag_orbit)

    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    simlog = synergia.utils.parallel_utils.Logger(0,
                synergia.utils.parallel_utils.LoggerV.ERROR)
    propagator.propagate(sim, simlog, 1)

    del propagator
    del simlog
    del stepper
    del diag_orbit
    del diag_bump
    del lp
    del bunch
    del sim

    h5_bump = h5py.File("bump.h5", 'r')
    h5_orbit = h5py.File("fulllattice.h5", 'r')

    bump_trks = h5_bump.get('track_coords')[()]
    orbit_trks = h5_orbit.get('track_coords')[()]
    h5_bump.close()
    h5_orbit.close()

    # test bump position
    assert target_x == pytest.approx(bump_trks[0, 0, 0])
    assert target_y == pytest.approx(bump_trks[0, 0, 2])

    # test particle reaches 0 at end of lattice
    assert orbit_trks[-1, 0, 0] == pytest.approx(0.0, abs=1.0e-8)
    assert orbit_trks[-1, 0, 1] == pytest.approx(0.0, abs=1.0e-8)
    assert orbit_trks[-1, 0, 2] == pytest.approx(0.0, abs=5.0e-7)
    assert orbit_trks[-1, 0, 3] == pytest.approx(0.0, abs=5.0e-8)

def main():
    test_three_bump()
    pass

if __name__ == "__main__":
    main()
