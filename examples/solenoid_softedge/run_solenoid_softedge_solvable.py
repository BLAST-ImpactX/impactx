#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

import amrex.space3d as amr
from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 250 MeV proton beam with an initial
# horizontal rms emittance of 1 um and an
# initial vertical rms emittance of 2 um
kin_energy_MeV = 250.0  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge

#   reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = 1.0 / 938.27208816 / 1e6  # electron charge/mass in e / eV

#   particle bunch
beam = sim.beam

dx = [1, 0, 0, 0, 0, 0]
dpx = [0, 1, 0, 0, 0, 0]
dy = [0, 0, 1, 0, 0, 0]
dpy = [0, 0, 0, 1, 0, 0]
dt = [0, 0, 0, 0, 1, 0]
dpt = [0, 0, 0, 0, 0, 1]

if amr.ParallelDescriptor.IOProcessor():
    beam.add_n_particles(dx, dy, dt, dpx, dpy, dpt, qm_eev, bunch_charge_C)

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

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

sim.lattice.extend(
    [
        monitor,
        sol,
        monitor,
    ]
)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
