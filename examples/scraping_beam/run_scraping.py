#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, Map6x6, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.max_level = 1
sim.n_cell = [16, 16, 20]
sim.blocking_factor_x = [16]
sim.blocking_factor_y = [16]
sim.blocking_factor_z = [4]

sim.space_charge = False
sim.poisson_solver = "fft"
sim.dynamic_size = True
sim.prob_relative = [1.2, 1.1]

# beam diagnostics
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 250  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 100000

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#   problem parameters
beam_radius = 2.0e-3
aperture_radius = 3.5e-3
correlation_k = 0.5
drift_distance = 6.0

#   particle bunch
distr = distribution.Kurth4D(
    lambdaX=beam_radius / 2.0,
    lambdaY=beam_radius / 2.0,
    lambdaT=beam_radius / 2.0,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# initialize the linear map
Iden = Map6x6.identity()
Rmat = Iden
Rmat[2, 1] = correlation_k
Rmat[4, 3] = correlation_k

# elements
drift1 = elements.Drift(
    name="d1",
    ds=drift_distance,
    aperture_x=aperture_radius,
    aperture_y=aperture_radius,
    nslice=40,
)
map1 = elements.LinearMap(R=Rmat)

# design the accelerator lattice
sim.lattice.extend([monitor, map1, drift1, monitor])

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
