#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import pandas as pd

import amrex.space3d as amr
from impactx import Config, ImpactX, elements

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
ref = sim.particle_container().push_ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

qm_eev = 1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV

pc = sim.particle_container()

if amr.ParallelDescriptor.IOProcessor():
    df_initial = pd.read_csv("./initial_coords.csv", sep=" ")
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
    if not Config.have_gpu:  # initialize using cpu-based PODVectors
        dx_podv = amr.PODVector_real_std()
        dy_podv = amr.PODVector_real_std()
        dt_podv = amr.PODVector_real_std()
        dpx_podv = amr.PODVector_real_std()
        dpy_podv = amr.PODVector_real_std()
        dpt_podv = amr.PODVector_real_std()
        dw_podv = amr.PODVector_real_std()
        dsx_podv = amr.PODVector_real_std()
        dsy_podv = amr.PODVector_real_std()
        dsz_podv = amr.PODVector_real_std()

    else:  # initialize on device using arena/gpu-based PODVectors
        dx_podv = amr.PODVector_real_arena()
        dy_podv = amr.PODVector_real_arena()
        dt_podv = amr.PODVector_real_arena()
        dpx_podv = amr.PODVector_real_arena()
        dpy_podv = amr.PODVector_real_arena()
        dpt_podv = amr.PODVector_real_arena()
        dw_podv = amr.PODVector_real_arena()
        dsx_podv = amr.PODVector_real_arena()
        dsy_podv = amr.PODVector_real_arena()
        dsz_podv = amr.PODVector_real_arena()

    for p_dx in dx:
        dx_podv.push_back(p_dx)
    for p_dy in dy:
        dy_podv.push_back(p_dy)
    for p_dt in dt:
        dt_podv.push_back(p_dt)
    for p_dpx in dpx:
        dpx_podv.push_back(p_dpx)
    for p_dpy in dpy:
        dpy_podv.push_back(p_dpy)
    for p_dpt in dpt:
        dpt_podv.push_back(p_dpt)
    for p_dw in dw:
        dw_podv.push_back(p_dw)
    for p_dsx in dsx:
        dsx_podv.push_back(p_dsx)
    for p_dsy in dsy:
        dsy_podv.push_back(p_dsy)
    for p_dsz in dsz:
        dsz_podv.push_back(p_dsz)

    pc.add_n_particles(
        dx_podv,
        dy_podv,
        dt_podv,
        dpx_podv,
        dpy_podv,
        dpt_podv,
        qm_eev,
        bunch_charge_C,
        None,
        dsx_podv,
        dsy_podv,
        dsz_podv,
    )

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5", period_sample_intervals=10)

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element
order = 4  # order of symplectic integration
nmap = 10  # number of steps for symplectic integration
fodo = [
    monitor,
    elements.ExactDrift(ds=0.25, nslice=ns),
    elements.ExactQuad(ds=1.0, k=1.0, int_order=order, mapsteps=nmap, nslice=ns),
    elements.ExactDrift(ds=0.5, nslice=ns),
    elements.ExactQuad(ds=1.0, k=-1.0, int_order=order, mapsteps=nmap, nslice=ns),
    elements.ExactDrift(ds=0.25, nslice=ns),
    monitor,
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
