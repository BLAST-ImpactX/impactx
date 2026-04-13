#!/usr/bin/python

# this script for Synergia2 and CHEF

import sys, os
import numpy as np

import synergia
import beamline

# create lattice with bend
lattice = synergia.lattice.Lattice("foo")

bend = synergia.lattice.Lattice_element("sbend", "bend")
# 90 degree bend with radius 1
bend.set_double_attribute("l", np.pi / 2)
bend.set_double_attribute("angle", np.pi / 2)

zdrift = synergia.lattice.Lattice_element("drift", "drift")
zdrift.set_double_attribute("l", 0)
lattice.append(zdrift)
lattice.append(bend)
lattice.append(zdrift)

mp = synergia.foundation.pconstants.mp

refpart = synergia.foundation.Reference_particle(1, mp, mp * 1.25)


lattice.set_reference_particle(refpart)

stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)

lattice_simulator = stepper.get_lattice_simulator()

chef_beamline = lattice_simulator.get_chef_lattice().get_beamline()

energy = chef_beamline.Energy()

print("energy: ", energy)

new_energy = energy + 0.0001
print("new energy: ", new_energy)

new_momentum = np.sqrt(new_energy**2 - mp**2)
old_momentum = np.sqrt(energy**2 - mp**2)

print("old_momentum: ", old_momentum)
print("new_momentum: ", new_momentum)
dpop = (new_momentum - old_momentum) / old_momentum

print("dpop: ", dpop)

p = beamline.Proton(energy)
p.set_ndp(dpop)

chef_beamline.propagate(p)

print("new cdt: ", p.get_cdt())
