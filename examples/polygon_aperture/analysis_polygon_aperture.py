#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#


import numpy as np
import openpmd_api as io

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = len(initial.momentum_x)
assert num_particles == len(initial)
assert num_particles != len(final)

print("Initial Beam: ", len(initial), " particles")
print("Final Beam: ", len(final), " particles")

# Make sure no particles are outside of the aperture in the final particle set
abs_x_final = abs(final["position_x"]).to_numpy()
abs_y_final = abs(final["position_y"]).to_numpy()

N = abs_x_final.shape[0]
insides = np.zeros(N, dtype=np.bool)
for i in range(N):
    insides[i] = ((abs_y_final[i] < 0.5e-3) and (abs_x_final[i] < 1.5e-3)) or (
        (abs_x_final[i] < 0.5e-3) and (abs_y_final[i] < 1.5e-3)
    )


outsides = ~insides
ninside = insides.sum()
noutside = outsides.sum()

assert ninside == len(final)
assert noutside == 0
