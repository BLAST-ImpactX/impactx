#!/usr/bin/env python

import sys

sys.path.append("..")
sys.path.append("../../foundation")
sys.path.append("../../lattice")
sys.path.append("../../simulation")
sys.path.append("../../bunch")
sys.path.append("/home/amundson/work/synergia2-devel_1_0/install/lib")

from nose.tools import *
from mad8_reader import Mad8_reader
from simulation import Lattice_simulator

from matching import generate_matched_bunch, generate_matched_bunch_transverse

# def test_3d():
#    lattice = Mad8_reader().get_lattice("fobrdobo", "foborodobo_s.lat")
#    lattice_simulator = Lattice_simulator(lattice, 1)
#    stdx = 0.001
#    stdy = 0.002
#    stdz = 0.01
#    num_real_particles = 1.0e12
#    num_macro_particles = 50
#    bunch = generate_matched_bunch(lattice_simulator, stdx, stdy, stdz,
#                           num_real_particles, num_macro_particles)


def test_transverse():
    lattice = Mad8_reader().get_lattice("fodo", "fodo.lat")
    #    lattice = Mad8_reader().get_lattice("fobrdobo", "foborodobo_s.lat")
    lattice_simulator = Lattice_simulator(lattice, 1)
    emit = 1e-6
    stdz = 0.01
    dpop = 1e-4

    num_real_particles = 1.0e12
    num_macro_particles = 50
    bunch = generate_matched_bunch_transverse(
        lattice_simulator,
        emit,
        emit,
        stdz,
        dpop,
        num_real_particles,
        num_macro_particles,
    )


#    assert(False)
