#!/usr/bin/env python

import sys

sys.path.append("../../..")
import local_paths

from synergia.optics import linear_one_turn_map
from synergia.lattice import Mad8_reader
from synergia.simulation import (
    Dummy_collective_operator,
    Lattice_simulator,
    Split_operator_stepper,
)


def test_linear_one_turn_map():
    num_steps = 1
    map_order = 2
    lattice = Mad8_reader().get_lattice("fodo", "lattices/fodo.lat")
    space_charge = Dummy_collective_operator("dummy")
    lattice_simulator = Lattice_simulator(lattice, map_order)
    stepper = Split_operator_stepper(lattice_simulator, space_charge, num_steps)
    map = linear_one_turn_map(lattice_simulator)
