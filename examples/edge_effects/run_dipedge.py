#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import pandas as pd

from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load initial beam parameters
kin_energy_MeV = 0.8e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = 1.0 / 938.27208816 / 1e6  # electron charge/mass in e / eV

beam = sim.beam

df_initial = pd.read_csv("./initial_coords.csv", sep=" ")
dx = df_initial["x"].to_numpy()
dpx = df_initial["px"].to_numpy()
dy = df_initial["y"].to_numpy()
dpy = df_initial["py"].to_numpy()
dt = df_initial["t"].to_numpy()
dpt = df_initial["pt"].to_numpy()
# dg = df_initial["gap"].to_numpy()
# print(dt)
print(dpt)
# print(dg)

beam.add_n_particles(dx, dy, dt, dpx, dpy, dpt, qm_eev, bunch_charge_C)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element
edge_angle = math.pi / 8.0
dipedge0 = elements.DipEdge(
    name="dipedge0", psi=edge_angle, rc=10.0, g=1.0e-3, model="linear", location="entry"
)
dipedge1 = elements.DipEdge(
    name="dipedge1",
    psi=edge_angle,
    rc=10.0,
    g=1.0e-3,
    model="nonlinear",
    location="entry",
)

line = [
    monitor,
    #    dipedge0,
    dipedge1,
    monitor,
]
# assign a segment
sim.lattice.extend(line)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
