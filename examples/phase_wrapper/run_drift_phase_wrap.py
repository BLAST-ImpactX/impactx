#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = True
# sim.particle_bc = "open"
# sim.particle_bc = "periodic"
# sim.particle_bc = "absorbing"
sim.particle_bc = "reflecting"

# domain decomposition & space charge mesh
sim.init_grids()

# load bunch parameters
kin_energy_MeV = 10.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 100000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

#   set bunch bucket length
sim.particle_container().set_bucket_length(0.23)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=3.9984884770e-5,
    lambdaY=3.9984884770e-5,
    lambdaT=5.0e-2,
    lambdaPx=2.6623538760e-5,
    lambdaPy=2.6623538760e-5,
    lambdaPt=0.0e-2,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 10  # number of slices per ds in the element
line = [
    monitor,
    elements.Drift(name="drift1", ds=0.5, nslice=ns),
    elements.ShortRF(name="shortrf1", V=10.0, freq=1.3e9, phase=-60.0),
    elements.Drift(name="drift1", ds=0.5, nslice=ns),
    monitor,
]
# assign the lattice
sim.lattice.extend(line)

# number of periods
sim.periods = 10

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
