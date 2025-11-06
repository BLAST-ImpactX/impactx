#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements, twiss

sim = ImpactX()

# set numerical parameters and IO control
sim.max_level = 0
sim.n_cell = [32, 32, 1]
sim.blocking_factor_x = [16]
sim.blocking_factor_y = [16]
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

# reference energy
kin_energy_MeV = 6.7  # reference energy

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.2720894).set_kin_energy_MeV(kin_energy_MeV)

#  beam current in A
beam_current_A = 0.5
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   particle bunch
distr = distribution.KVdist(
    **twiss(
        beta_x=0.737881,
        beta_y=0.737881,
        beta_t=0.5,
        emitt_x=1.0e-6,
        emitt_y=1.0e-6,
        emitt_t=1.0e-12,
        alpha_x=2.4685083,
        alpha_y=-2.4685083,
        alpha_t=0.0,
    )
)
sim.add_particles(beam_current_A, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 100  # number of slices per ds in the element
fodo = [
    monitor,
    elements.Drift(name="drift1", ds=7.44e-2, nslice=ns),
    elements.Quad(name="quad1", ds=6.10e-2, k=-103.12574100336, nslice=ns),
    elements.Drift(name="drift2", ds=14.88e-2, nslice=ns),
    elements.Quad(name="quad2", ds=6.10e-2, k=103.12574100336, nslice=ns),
    elements.Drift(name="drift3", ds=7.44e-2, nslice=ns),
    monitor,
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
