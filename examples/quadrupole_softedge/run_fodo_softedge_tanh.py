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

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=3.9984884770e-5,
    lambdaY=3.9984884770e-5,
    lambdaT=1.0e-3,
    lambdaPx=2.6623538760e-5,
    lambdaPy=2.6623538760e-5,
    lambdaPt=2.0e-3,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 10  # number of slices per ds in the element

# specify the on-axis field profile
zmin = -0.75  # lower value of on-axis longitudinal coordinate (in meters)
zmax = 0.75  # upper value of on-axis longitudinal coordinate (in meters)
nz = 401  # number of longitudinal sampling points to be used
g = 0.07  # gap parameter (in meters)
L_hardedge = 1.0  # length in the hard-edge limit
zdata = np.linspace(zmin, zmax, nz)
bdata = (
    1.0
    / 2.0
    * (
        np.tanh((zdata + L_hardedge / 2.0) / g)
        - np.tanh((zdata - L_hardedge / 2.0) / g)
    )
)

# hard-edge lattice
quad1_he = elements.Quad(name="quad1", ds=1.0, k=1.0, nslice=ns)
quad2_he = elements.Quad(name="quad2", ds=1.0, k=-1.0, nslice=ns)

# design the accelerator lattice
quad1 = elements.SoftQuadrupole(
    name="quad1",
    ds=1.5,
    gscale=1.0,
    z=zdata,
    gradient_on_axis=bdata,
    ncoef=35,
    mapsteps=200,
    nslice=ns,
)

# design the accelerator lattice
quad2 = elements.SoftQuadrupole(
    name="quad1",
    ds=1.5,
    gscale=-1.0,
    z=zdata,
    gradient_on_axis=bdata,
    ncoef=35,
    mapsteps=200,
    nslice=ns,
)

# hard-edge lattice
drift1_he = elements.Drift(name="drift1", ds=0.25, nslice=ns)
drift2_he = elements.Drift(name="drift2", ds=0.5, nslice=ns)

# soft-edge lattice
drift1 = elements.Drift(name="drift1", ds=0.0, nslice=ns)
drift2 = elements.Drift(name="drift2", ds=0.0, nslice=ns)

# assign a fodo segment
sim.lattice.extend([monitor, drift1, quad1, drift2, quad2, drift1, monitor])
# sim.lattice.extend([monitor, drift1_he, quad1_he, drift2_he, quad2_he, drift1_he, monitor])

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
