#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
beam_initial = series.iterations[1].particles["beam"]
initial_sort = beam_initial.to_df().set_index("id")
beam_final = series.iterations[last_step].particles["beam"]
final_sort = beam_final.to_df().set_index("id")

# Sextupole parameters
K_normal = 10.0
K_skew = -2.0

# Initial particle data

xi = initial_sort["position_x"]
pxi = initial_sort["momentum_x"]
yi = initial_sort["position_y"]
pyi = initial_sort["momentum_y"]
ti = initial_sort["position_t"]
pti = initial_sort["momentum_t"]

# Predicted momentum kick, from MAD-X user guide (equation 1.15)

px_predicted = K_normal / 2.0 * (yi**2 - xi**2) + K_skew * xi * yi
py_predicted = K_skew / 2.0 * (xi**2 - yi**2) + K_normal * xi * yi

# Maximum momentum kick
px_max = px_predicted.abs().max()
py_max = py_predicted.abs().max()

# Final particle data

xf = final_sort["position_x"]
pxf = final_sort["momentum_x"]
yf = final_sort["position_y"]
pyf = final_sort["momentum_y"]

# Difference between value and prediction

dpx_max = (pxf - px_predicted).abs().max()
dpy_max = (pyf - py_predicted).abs().max()

print()
print("Difference between predicted and computed final momentum, absolute max:")
print("dpx_max", dpx_max)
print("dpy_max", dpy_max)

# Test maximum error:
atol = 5.1e-12
print(f"  tol={atol}")

print()
print("Difference between predicted and computed final momentum (max), relative:")
print("dpx_max/px_max", dpx_max / px_max)
print("dpy_max/py_max", dpy_max / py_max)

# Test maximum error:
atol = 5.1e-12
print(f"  tol={atol}")

assert np.allclose(
    [dpx_max / px_max, dpy_max / py_max],
    [0.0, 0.0],
    atol=atol,
)
