#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import pandas as pd

import amrex.space3d as amr
from impactx import Config, ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.max_level = 0
sim.n_cell = [32, 32, 1]
sim.blocking_factor_x = [32]
sim.blocking_factor_y = [32]
sim.blocking_factor_z = [1]

sim.particle_shape = 2  # B-spline order
sim.space_charge = "2D"
sim.poisson_solver = "fft"
sim.dynamic_size = True
sim.prob_relative = [1.1]

# beam diagnostics
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 250  # reference energy
beam_current_A = 0.15  # beam current
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = 1.0 / 938.27208816 / 1e6  # electron charge/mass in e / eV

#   particle bunch
distr = distribution.KVdist(
    lambdaX=5.0e-4,
    lambdaY=5.0e-4,
    lambdaT=1.0e-3,
    lambdaPx=0.0,
    lambdaPy=0.0,
    lambdaPt=0.0,
)
sim.add_particles(beam_current_A, distr, npart)

pc = sim.particle_container()

#  add test particles
df_initial = pd.read_csv("./initial_coords.csv", sep=" ")
dx = df_initial["x"].to_numpy()
dpx = df_initial["px"].to_numpy()
dy = df_initial["y"].to_numpy()
dpy = df_initial["py"].to_numpy()
dt = df_initial["t"].to_numpy()
dpt = df_initial["pt"].to_numpy()
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

pc.add_n_particles(dx_podv, dy_podv, dt_podv, dpx_podv, dpy_podv, dpt_podv, qm_eev, 0.0)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
doubling_distance = 10.612823669911099

sim.lattice.extend(
    [monitor, elements.Drift(name="d1", ds=doubling_distance, nslice=100), monitor]
)

sarr = []
test_data = []
mm_scale = 1.0e3


def hook_before_slice(sim):
    s = sim.beam.ref.s
    sarr.append(s)
    beam = sim.beam.to_df()
    # Filter on particle weight (collect test particles only)
    for row in beam[beam["weighting"] == 0.0].itertuples():
        # collect test particle data
        test_data.append(
            [s, row.idcpu, row.position_x * mm_scale, row.position_y * mm_scale]
        )


sim.hook["before_slice"] = hook_before_slice

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()

df = pd.DataFrame(test_data, columns=["s", "id", "x", "y"])
sorted_df = df.sort_values(by="id")

n = len(sarr)
for i in range(0, len(df), n):
    subset = sorted_df.iloc[i : i + n]
    plt.scatter(subset["s"], subset["x"], s=5)

plt.xlabel("s [m]", fontsize=12)
plt.ylabel("x [mm]", fontsize=12)
plt.title("Test Particles: Horizontal Coordinates")
plt.show()
