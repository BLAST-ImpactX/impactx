#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
import pandas as pd

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = 7
assert num_particles == len(initial)
assert num_particles == len(final)

# load particle data
df_initial = pd.read_csv("./initial_coords.csv", sep=" ")
df_final = pd.read_csv("./final_coords.csv", sep=" ")

# initial coordinates
xi = df_initial["x"]
pxi = df_initial["px"]
yi = df_initial["y"]
pyi = df_initial["py"]
ti = df_initial["t"]
pti = df_initial["pt"]
sxi = df_initial["sx"]
syi = df_initial["sy"]
szi = df_initial["sz"]
norm_phase_space_vector = np.sqrt(xi**2 + pxi**2 + yi**2 + pyi**2 + ti**2 + pti**2)

# difference of final coordinates from linear values
diff_xf = df_final["x"] - final["position_x"]
diff_pxf = df_final["px"] - final["momentum_x"]
diff_yf = df_final["y"] - final["position_y"]
diff_pyf = df_final["py"] - final["momentum_y"]
diff_tf = df_final["t"] - final["position_t"]
diff_ptf = df_final["pt"] - final["momentum_t"]
diff_sxf = df_final["sx"] - final["spin_x"]
diff_syf = df_final["sy"] - final["spin_y"]
diff_szf = df_final["sz"] - final["spin_z"]
norm_spin_diff = np.sqrt(diff_sxf**2 + diff_syf**2 + diff_szf**2)

print("")
print(":")
print("Phase space vector norm:")
print(norm_phase_space_vector)

print(":")
print("Spin difference norm:")
print(norm_spin_diff)

m, b = np.polyfit(np.log(norm_phase_space_vector), np.log(norm_spin_diff), deg=1)
print(f"  slope={m}")

# Test for  quadratic scaling with initial phase space vector:
rtol = 2.0e-2
print(f"  rtol={rtol}")
assert np.allclose(
    [m],
    [
        2.0,
    ],
    rtol=rtol,
)
