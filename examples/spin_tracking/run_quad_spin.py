#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.spin = True
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# basic beam parameters
kin_energy_MeV = 100.0  # reference energy (kinetic)
bunch_charge_C = 25.0e-12  # used with space charge
npart = 100000  # number of macro particles

# set reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=0.0003,
    lambdaY=0.0003,
    lambdaT=1.0e-6,
    lambdaPx=0.0002,
    lambdaPy=0.0002,
    lambdaPt=1.0e-6,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
spin = distribution.SpinvMF(
    0.4,
    0.9,
    0.1,
)

sim.add_particles(bunch_charge_C, distr, npart, spin)

# design the accelerator lattice
ns = 1  # number of slices per ds in the element

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# quad focusing strength (> 0)
k_value = 100.0
print("k_value")
print(k_value)

# length for this test should be one period
ds_value = 2.0 * np.pi / np.sqrt(k_value)

quad1 = elements.Quad(name="quad1", ds=ds_value, k=k_value, nslice=ns)

# set the lattice
sim.lattice.append(monitor)
sim.lattice.append(quad1)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
