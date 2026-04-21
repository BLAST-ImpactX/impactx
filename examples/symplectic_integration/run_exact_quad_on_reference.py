#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Nikita Kuklev, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


import amrex.space3d as amr
from impactx import Config, ImpactX, elements

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

dx = [0.0]
dpx = [0.0]
dy = [0.0]
dpy = [0.0]
dt = [0.0]
dpt = [0.0]
if not Config.have_gpu:  # initialize using cpu-based PODVectors
    dx_podv = amr.PODVector_real_std()
    dy_podv = amr.PODVector_real_std()
    dt_podv = amr.PODVector_real_std()
    dpx_podv = amr.PODVector_real_std()
    dpy_podv = amr.PODVector_real_std()
    dpt_podv = amr.PODVector_real_std()
else:  # initialize on device using arena/gpu-based PODVectors
    dx_podv = amr.PODVector_real_arena()
    dy_podv = amr.PODVector_real_arena()
    dt_podv = amr.PODVector_real_arena()
    dpx_podv = amr.PODVector_real_arena()
    dpy_podv = amr.PODVector_real_arena()
    dpt_podv = amr.PODVector_real_arena()

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

beam.add_n_particles(
    dx_podv, dy_podv, dt_podv, dpx_podv, dpy_podv, dpt_podv, qm_eev, bunch_charge_C
)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element
line = [
    monitor,
    elements.Quad(
        name="q",
        ds=0.5,
        k=1.2,
        nslice=ns,
    ),
    elements.ExactQuad(
        name="eq",
        ds=0.5,
        k=1.2,
        nslice=ns,
    ),
    elements.ExactMultipole(
        name="em",
        ds=0.5,
        k_normal=[0.0, 1.2],
        k_skew=[0.0, 0.0],
        unit=0,
        int_order=4,
        mapsteps=5,
        nslice=ns,
    ),
    monitor,
]

# assign a fodo segment
sim.lattice.extend(line)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
