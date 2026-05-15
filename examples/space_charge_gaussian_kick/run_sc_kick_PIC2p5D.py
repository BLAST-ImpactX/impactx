#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Ji Qiang
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = "2p5D"
sim.poisson_solver = "fft"
sim.space_charge_num_longitudinal_bins = 129
sim.space_charge_apply_longitudinal_kick = True
sim.slice_step_diagnostics = True

# set numerical parameters and IO control
sim.max_level = 0
sim.n_cell = [64, 64, 1]
sim.blocking_factor_x = [16]
sim.blocking_factor_y = [16]
sim.blocking_factor_z = [1]
sim.prob_relative = [1.1]

# beam diagnostics
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()


# domain decomposition & space charge mesh
sim.init_grids()

# beam energy ad bunch charge
kin_energy_MeV = 100  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 100000  # number of macro particles

# intialize reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

# initialize particle bunch
distr = distribution.Gaussian(
    lambdaX=4.0e-5,
    lambdaY=4.0e-5,
    lambdaT=1.0e-3,
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
L = 1.0  # drift length
sim.lattice.extend([monitor, elements.Drift(name="d1", ds=L, nslice=ns), monitor])

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
