#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from scipy.constants import c, eV, m_p

from impactx import ImpactX, distribution, elements

mp_mev = 1.0e-6 * m_p * c**2 / eV
total_Booster_charge = 6.7e12  # PIP-II full Booster
active_buckets = 81  # 81 out of 84 buckets full

turns = 1

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 800 MeV proton beam

kin_energy_MeV = 800.0  # reference energy 800 MeV
bunch_charge_C = eV * total_Booster_charge / active_buckets  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

# this is not the distribution that will be used
#   particle bunch
distr = distribution.Waterbag(
    lambdaX=3.9984884770e-5,
    lambdaY=3.9984884770e-5,
    lambdaT=1.0e-3,
    lambdaPx=2.6623538760e-5,
    lambdaPy=2.6623538760e-5,
    lambdaPt=2.0e-3,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# Read the Booster lattice
with open("booster_impactx_lattice.txt", "r") as F:
    lattice_txt = F.read()
booster = eval(lattice_txt)
sim.lattice.extend(booster)
sim.lattice.append(monitor)

# run simulation
sim.periods = turns
sim.track_particles()

# clean shutdown
sim.finalize()
