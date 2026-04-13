#!/usr/bin/env python

import numpy
from synergia.lattice import reference_particle_to_chef_jet_particle, \
    chef_unit_conversion, get_chef_index

import basic_toolkit
import mxyzptlk
from beamline import JetProton, JetParticle

def _convert_linear_map(chef_map, reference_particle):
    u = chef_unit_conversion(reference_particle)
    map = numpy.zeros((6, 6), 'd')
    for row in range(0, 6):
        for column in range(0, 6):
            chef_row = get_chef_index(row)
            chef_column = get_chef_index(column)
            map[row, column] = chef_map.get(chef_row, chef_column) * \
                              u[row] / u[column]
    return map

def linear_one_turn_map(lattice_simulator):
    reference_particle = lattice_simulator.get_lattice().get_reference_particle()
    map_order = lattice_simulator.get_map_order()
    JetParticle.createStandardEnvironments(map_order)
    jet_particle = JetProton(reference_particle.get_four_momentum().get_total_energy())
   # jet_particle = reference_particle_to_chef_jet_particle(reference_particle,
   #                                                        map_order)
    lattice_simulator.get_chef_lattice().get_beamline().propagateJetParticle(jet_particle)    
    return _convert_linear_map(jet_particle.State().jacobian(), reference_particle)


