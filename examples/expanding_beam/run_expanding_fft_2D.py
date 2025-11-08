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
sim.max_level = 0
sim.n_cell = [32, 32, 1]
sim.blocking_factor_x = [32]
sim.blocking_factor_y = [32]
sim.blocking_factor_z = [1]

sim.particle_shape = 2  # B-spline order
sim.space_charge = "2D"
sim.poisson_solver = "fft"
sim.dynamic_size = True
sim.prob_relative = [1.1]

# beam diagnostics
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 250  # reference energy
beam_current_A = 0.15  # beam current
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.2720894).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.KVdist(
    lambdaX=5.0e-4,
    lambdaY=5.0e-4,
    lambdaT=1.0e-3,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
)
sim.add_particles(beam_current_A, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
doubling_distance = 10.612823669911099

sim.lattice.extend(
    [monitor, elements.Drift(name="d1", ds=doubling_distance, nslice=100), monitor]
)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
