#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, Map3x6, Map6x6, Vector3, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# reference kinetic energy
kin_energy_MeV = 250.0  # reference energy

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

# specify the on-axis field profile
zmin = -1.0  # lower value of on-axis longitudinal coordinate (in meters)
zmax = 1.0  # upper value of on-axis longitudinal coordinate (in meters)
nz = 401  # number of longitudinal sampling points to be used
g = 0.1  # gap parameter (in meters)
zdata = np.linspace(zmin, zmax, nz)
bdata = 1.0 / (1.0 + (zdata / g) ** 2)

# design the accelerator lattice
sol = elements.SoftSolenoid(
    name="sol1",
    ds=2.0,
    bscale=-1.0,
    z=zdata,
    field_on_axis=bdata,
    ncoef=35,
    mapsteps=200,
    nslice=1,
)

# lattice
sim.lattice.extend(
    [
        sol,
    ]
)

# run simulation
sim.track_reference(ref)

# return linear map
Rmat = Map6x6()
vmat = Vector3()
Amat = Map3x6()

Rmat = ref.map
vmat = ref.spin_rotation_vector
Amat = ref.spin_coupling

print()
print("Linear map:")
for i in range(1, 7):
    for j in range(1, 7):
        print(i, j, Rmat[i, j])

print()
print("Reference spin rotation vector:")
for i in range(1, 4):
    print(i, vmat[i])

print()
print("Spin-orbit coupling matrix:")
for i in range(1, 4):
    for j in range(1, 7):
        print(i, j, Amat[i, j])

# clean shutdown
sim.finalize()
