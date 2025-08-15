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
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load initial beam
kin_energy_MeV = 1.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=1.0e-3,
    lambdaY=1.0e-3,
    lambdaT=0.3,
    lambdaPx=2.0e-4,
    lambdaPy=2.0e-4,
    lambdaPt=2.0e-5,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# problem parameters

rigidity = 3.337345025729098
g_gradient = 100.0
k_gradient = g_gradient / rigidity
length = 0.25
angle_deg = 45.0
angle_rad = np.pi * angle_deg / 180.0
field = rigidity * angle_rad / length

# design the accelerator lattice)
ns = 20  # number of slices per ds in the element

quad1 = elements.ChrQuad(name="quad1", ds=length, k=k_gradient, unit=0, nslice=ns)
quad2inv = elements.ChrQuad(
    name="quad2inv", ds=-length, k=-g_gradient, unit=1, nslice=ns
)

quad3 = elements.ChrQuad(name="quad3", ds=length, k=-k_gradient, unit=0, nslice=ns)
quad4inv = elements.ChrQuad(
    name="quad4inv", ds=-length, k=g_gradient, unit=1, nslice=ns
)

bend1 = elements.ExactSbend(name="bend1", ds=length, phi=angle_deg, B=field, nslice=ns)
bend2inv = elements.ExactSbend(name="bend2inv", ds=-length, phi=angle_deg, nslice=ns)

bend3 = elements.ExactSbend(name="bend3", ds=length, phi=angle_deg, B=-field, nslice=ns)
bend4inv = elements.ExactSbend(name="bend4inv", ds=-length, phi=-angle_deg, nslice=ns)

line = [
    monitor,
    quad1,
    quad2inv,
    quad3,
    quad4inv,
    bend1,
    bend2inv,
    bend3,
    bend4inv,
    monitor,
]

sim.lattice.extend(line)

# number of periods through the lattice
sim.periods = 1

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
