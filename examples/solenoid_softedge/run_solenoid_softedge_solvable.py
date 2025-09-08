#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import amrex.space3d as amr
from impactx import Config, ImpactX, elements

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
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)
qm_eev = 1.0 / 938.27208816 / 1e6  # electron charge/mass in e / eV

#   particle bunch
pc = sim.particle_container()

dx = [1, 0, 0, 0, 0, 0]
dpx = [0, 1, 0, 0, 0, 0]
dy = [0, 0, 1, 0, 0, 0]
dpy = [0, 0, 0, 1, 0, 0]
dt = [0, 0, 0, 0, 1, 0]
dpt = [0, 0, 0, 0, 0, 1]

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

if amr.ParallelDescriptor.IOProcessor():
    pc.add_n_particles(
        dx_podv, dy_podv, dt_podv, dpx_podv, dpy_podv, dpt_podv, qm_eev, bunch_charge_C
    )

# design the accelerator lattice
sol = elements.SoftSolenoid(
    name="sol1",
    ds=2.0,
    bscale=-1.0,
    cos_coefficients=[
        0.294225534860747,
        0.2317701877038734,
        0.166794824926630,
        0.122811705237718,
        0.089180156215475,
        0.065459443762749,
        0.047593623793664,
        0.034919898089417,
        0.0253867332415233,
        0.018635438261000,
        0.013536778750496,
        0.009948507321628,
        0.0072152939652668,
        0.005313407598416,
        0.003843724102950,
        0.002839750408934,
        0.0020458821109883,
        0.001519293408286,
        0.001087477417341,
        0.000814190355962,
        0.000576771050997,
    ],
    sin_coefficients=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
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
