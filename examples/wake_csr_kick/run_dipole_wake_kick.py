#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import eV

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.short_range_wakefields = True
sim.wake_bins = 150
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# parameters specific to the Booster
total_Booster_charge = 6.7e12  # PIP-II full Booster
active_buckets = 81  # 81 out of 84 buckets full

# load a 800 MeV proton beam
kin_energy_MeV = 800.0  # reference energy 800 MeV
bunch_charge_C = eV * total_Booster_charge / active_buckets  # used with space charge
npart = 1000000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.KVdist(
    lambdaX=0.0,
    lambdaY=0.0,
    lambdaT=1.0e-3,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice

# access RF cavity on-axis field data
data_in = np.loadtxt("onaxis_data.in")
z = data_in[:, 0]
ez_onaxis = data_in[:, 1]
ncoef = 25

# specify the RF cavity short-range structure wake data using an expression
zmin = (
    0.0  # lower value of longitudinal separation used in the wake function (in meters)
)
zmax = (
    0.01  # upper value of longitudinal separation used in the wake function (in meters)
)
nz = 5000  # number of longitudinal separation sampling points to be used
K0 = 4.30355e13  # first fitted parameter
s0 = 0.00175  # second fitted parameter (in meters)
z_wake = np.linspace(zmin, zmax, nz)
wake_z_data = K0 * np.exp(-np.sqrt(z_wake / s0))
# This analytical functional form is fitted to the tabulated single-particle wake for a TESLA, 9-cell 1.3 GHz cavity as it appears in (https://github.com/impact-lbl/IMPACT-Z/tree/master/examples/Example3).

#   Focusing combined-function dipole element
cfb = elements.ExactCFbend(
    ds=2.889612,
    int_order=6,
    k_normal=[0.024481550532148122, 0.05410921561, 0.00546850017312],
    k_skew=[0.0, 0.0, 0.0],   
    z_data=z_wake,
    wake_z_data=wake_z_data,
    mapsteps=6,
    unit=0,
    nslice=1,
)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

sim.lattice.extend(
    [
        monitor,
        cfb,
        monitor,
    ]
)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
