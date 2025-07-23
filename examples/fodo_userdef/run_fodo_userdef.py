#!/usr/bin/env python3
#
# Copyright 2022-2024 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

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
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    **twiss(
        beta_x=2.8216194100262637,
        beta_y=2.8216194100262637,
        beta_t=0.5,
        emitt_x=2e-09,
        emitt_y=2e-09,
        emitt_t=2e-06,
        alpha_x=-1.5905003499999992,
        alpha_y=1.5905003499999992,
        alpha_t=0.0,
    )
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# add a user-defined, linear element for the drifts
Iden = Map6x6.identity()
R1, R2 = Iden, Iden

ds1 = 0.25
R1[1, 2] = ds1
R1[3, 4] = ds1
R1[5, 6] = ds1 / 16.6464  # ds / (beta*gamma^2)
drift1 = elements.LinearMap(name="drift1", R=R1, ds=ds1)

ds2 = 0.5
R2[1, 2] = ds2
R2[3, 4] = ds2
R2[5, 6] = ds2 / 16.6464  # ds / (beta*gamma^2)
drift2 = elements.LinearMap(name="drift2", R=R2, ds=ds2)

# design the accelerator lattice)
ns = 25  # number of slices per ds in the element
fodo = [
    monitor,
    drift1,
    monitor,
    elements.Quad(name="quad1", ds=1.0, k=1.0, nslice=ns),
    monitor,
    drift2,
    monitor,
    elements.Quad(name="quad2", ds=1.0, k=-1.0, nslice=ns),
    monitor,
    drift1,
    monitor,
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
