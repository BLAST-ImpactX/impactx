#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.max_level = 0
sim.n_cell = [32, 32, 1]
sim.blocking_factor_x = [4]
sim.blocking_factor_y = [4]
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

# basic beam properties
R0 = 1.0e-3  # initial beam radius (m)
Ib = 0.15  # beam current (A)
kin_energy_MeV = 250  # reference energy
beam_current_A = Ib  # beam current
npart = 10000  # number of macro particles (outside tests, use 1e5 or more)

#   reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = 1.0 / 938.27208816 / 1e6  # electron charge/mass in e / eV

#   particle bunch
distr = distribution.KVdist(
    lambdaX=R0 / 2.0,
    lambdaY=R0 / 2.0,
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
pc.add_n_particles(dx, dy, dt, dpx, dpy, dpt, qm_eev, bunch_charge=0.0)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
betgam = ref.beta_gamma  # relativistic factor
IA = 3.1297388031196285e7  # alfven current for protons
Kpv = 2.0 * Ib / (IA * betgam**3)  # generalized beam perveance
double_constant = 1.516770632602484  # constant independent of beam parameters
doubling_distance = R0 * double_constant / np.sqrt(Kpv)
print("Doubling distance: ")
print(doubling_distance)

ns = 50  # recommended: for outside tests, using 100-200.

sim.lattice.extend(
    [monitor, elements.Drift(name="d1", ds=doubling_distance, nslice=ns), monitor]
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
