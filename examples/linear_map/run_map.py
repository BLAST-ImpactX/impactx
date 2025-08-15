#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

# from elements import LinearTransport
import numpy as np

from impactx import ImpactX, Map6x6, distribution, elements, twiss

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 45.6e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   target beta functions (m)
beta_star_x = 0.15
beta_star_y = 0.8e-3
beta_star_t = 9.210526315789473

#   particle bunch
distr = distribution.Waterbag(
    **twiss(
        beta_x=beta_star_x,
        beta_y=beta_star_y,
        beta_t=beta_star_t,
        emitt_x=0.27e-09,
        emitt_y=1.0e-12,
        emitt_t=1.33e-06,
        alpha_x=0.0,
        alpha_y=0.0,
        alpha_t=0.0,
    )
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# initialize the linear map
Iden = Map6x6.identity()
Rmat = Iden

# desired tunes
Qx = 0.139
Qy = 0.219
Qt = 0.0250

# desired phase advance
phi_x = 2.0 * np.pi * Qx
phi_y = 2.0 * np.pi * Qy
phi_t = 2.0 * np.pi * Qt

# matrix elements for the horizontal plane
Rmat[1, 1] = np.cos(phi_x)
Rmat[1, 2] = beta_star_x * np.sin(phi_x)
Rmat[2, 1] = -np.sin(phi_x) / beta_star_x
Rmat[2, 2] = np.cos(phi_x)
# matrix elements for the vertical plane
Rmat[3, 3] = np.cos(phi_y)
Rmat[3, 4] = beta_star_y * np.sin(phi_y)
Rmat[4, 3] = -np.sin(phi_y) / beta_star_y
Rmat[4, 4] = np.cos(phi_y)
# matrix elements for the longitudinal plane
Rmat[5, 5] = np.cos(phi_t)
Rmat[5, 6] = beta_star_t * np.sin(phi_t)
Rmat[6, 5] = -np.sin(phi_t) / beta_star_t
Rmat[6, 6] = np.cos(phi_t)

# design the accelerator lattice
map = [
    monitor,
    elements.LinearMap(R=Rmat),
]

sim.lattice.extend(map)

# number of periods through the lattice
sim.periods = 4

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
