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
kin_energy_MeV = 2000.0  # reference energy (kinetic)
mass_MeV = 0.510998950  # particle mass
bunch_charge_C = 25.0e-12  # used with space charge
gyromagnetic_anomaly = 0.00115965218062  # value for an electron
npart = 100000  # number of macro particles

# set reference particle
ref = sim.particle_container().push_ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=0.003,
    lambdaY=0.003,
    lambdaT=1.0e-6,
    lambdaPx=0.2,
    lambdaPy=0.2,
    lambdaPt=0.2,
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

# bend radius (> 0)
rc_value = 1.0
print("rc_value")
print(rc_value)

gamma = kin_energy_MeV / mass_MeV + 1
print("relativistic gamma")
print(gamma)

# length for this test should be one period
ds_value = 2.0 * np.pi / (gyromagnetic_anomaly * gamma)

bend1 = elements.Sbend(name="bend1", ds=ds_value, rc=rc_value, nslice=ns)

# set the lattice
sim.lattice.append(monitor)
sim.lattice.append(bend1)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
