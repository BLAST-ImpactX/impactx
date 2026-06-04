#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

import amrex.space3d as amr
from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.spin = True
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# beam parameters
kin_energy_MeV = 10.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

qm_eev = 1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV

beam = sim.beam

if amr.ParallelDescriptor.IOProcessor():
    df_initial = pd.read_csv("./initial_coords_sbend.csv", sep=" ")
    dx = df_initial["x"].to_numpy()
    dpx = df_initial["px"].to_numpy()
    dy = df_initial["y"].to_numpy()
    dpy = df_initial["py"].to_numpy()
    dt = df_initial["t"].to_numpy()
    dpt = df_initial["pt"].to_numpy()
    dw = df_initial["w"].to_numpy()
    dsx = df_initial["sx"].to_numpy()
    dsy = df_initial["sy"].to_numpy()
    dsz = df_initial["sz"].to_numpy()
    beam.add_n_particles(
        dx, dy, dt, dpx, dpy, dpt, qm_eev, bunch_charge_C, sx=dsx, sy=dsy, sz=dsz
    )

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5", period_sample_intervals=10)

# design the accelerator lattice
ds_value = 1.0
rc_value = 10.0
phi_value = 180.0 / np.pi * (ds_value / rc_value)
ns = 1

# to use the magnetic field instead
Brho = ref.rigidity_Tm
B_value = Brho / rc_value

bend1 = elements.ExactSbend(name="bend1", ds=ds_value, phi=phi_value, nslice=ns)
bend2 = elements.Sbend(name="bend2", ds=-ds_value, rc=rc_value, nslice=ns)
bend3 = elements.ExactSbend(
    name="bend3", ds=ds_value, phi=phi_value, B=B_value, nslice=ns
)

# set the lattice
sim.lattice.append(monitor)

# test of forward + reverse using rc = ds/phi
sim.lattice.append(bend1)
sim.lattice.append(bend2)

# test of forward + reverse using rc = Brho/B
sim.lattice.append(bend3)
sim.lattice.append(bend2)

sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
