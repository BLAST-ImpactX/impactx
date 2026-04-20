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
gyromagnetic_anomaly = 2.0

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
ref.set_gyromagnetic_anomaly(gyromagnetic_anomaly)

# specify the on-axis field profile
zmin = -1.0  # lower value of on-axis longitudinal coordinate (in meters)
zmax = 1.0  # upper value of on-axis longitudinal coordinate (in meters)
# nz = 401  # number of longitudinal sampling points to be used
nz = 801
g = 1.0  # gap parameter (in meters)
zdata = np.linspace(zmin, zmax, nz)
bdata = 1.0 / (1.0 + (zdata / g) ** 2)
bscale = 5.0 / 6.0

# design the accelerator lattice
sol = elements.SoftSolenoid(
    name="sol1",
    ds=2.0,
    bscale=bscale,
    z=zdata,
    field_on_axis=bdata,
    #    ncoef=35,
    #    mapsteps=800,
    ncoef=50,
    mapsteps=800,
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

# return spin map
Rmat = Map6x6()
vmat = Vector3()
Amat = Map3x6()
vpred = Vector3()
Apred = Map3x6()

Rmat = ref.map  # not used - included for illustration
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

gamma = ref.gamma
beta = ref.beta
vpred[3] = -5.0 * np.pi / 4.0

Apred[1, 1] = (
    (56354 + 8125 * np.pi - 5876 * np.sqrt(2.0) - 6676 * np.sqrt(6))
    * (gamma - 1.0)
    / (81120 * g)
)
Apred[1, 2] = (
    (2146 + 325 * np.pi + 52 * np.sqrt(2.0) - 500 * np.sqrt(6)) * (gamma - 1.0) / (3380)
)
Apred[1, 3] = (
    (6850 + 325 * np.pi - 6676 * np.sqrt(2.0) + 5876 * np.sqrt(6))
    * (gamma - 1.0)
    / (81120 * g)
)
Apred[1, 4] = (
    -(4146 + 325 * np.pi + 500 * np.sqrt(2.0) + 52 * np.sqrt(6))
    * (gamma - 1.0)
    / (3380)
)
Apred[3, 6] = -5 * np.pi / (4.0 * beta)
Apred[2, 1] = -Apred[1, 3]
Apred[2, 2] = -Apred[1, 4]
Apred[2, 3] = Apred[1, 1]
Apred[2, 4] = Apred[1, 2]

print()
print("Reference spin rotation vector, predicted:")
for i in range(1, 4):
    print(i, vpred[i])

print()
print("Spin-orbit coupling matrix, predicted:")
for i in range(1, 4):
    for j in range(1, 7):
        print(i, j, Apred[i, j])

# clean shutdown
sim.finalize()

# analysis
atol = 2.0e-6
rtol = 0.0
print(f"  atol={atol}")

assert np.allclose(vmat, vpred, rtol=rtol, atol=atol)

atol = 2.0e-4  # can this tolerance be improved?
rtol = 0.0
print(f"  atol={atol}")

assert np.allclose(Amat, Apred, rtol=rtol, atol=atol)
