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
sim.particle_bc = "periodic"
# sim.particle_bc = "absorbing"
# sim.particle_bc = "reflecting"

# domain decomposition & space charge mesh
sim.init_grids()

# load bunch parameters
kin_energy_MeV = 10.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   set bunch bucket length
sim.particle_container().set_bucket_length(0.2)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=5.0e-4,
    lambdaY=5.0e-4,
    lambdaT=0.2,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element
line = [
    monitor,
    monitor,
]
# assign the lattice
sim.lattice.extend(line)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
