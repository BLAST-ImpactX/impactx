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
initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = 6
assert num_particles == len(initial)
assert num_particles == len(final)

# initial data
xi = np.array([1, 0, 0, 0, 0, 0])
pxi = np.array([0, 1, 0, 0, 0, 0])
yi = np.array([0, 0, 1, 0, 0, 0])
pyi = np.array([0, 0, 0, 1, 0, 0])
ti = np.array([0, 0, 0, 0, 1, 0])
pti = np.array([0, 0, 0, 0, 0, 1])

# problem parameters
g = 0.1
L = 2.0
bscale = -1.0  # This coincides with the value set in the input file.
gamma_ref = beam_final.get_attribute("gamma_ref")

# derived parameters
lbda = L / (2.0 * g)
phi = 2.0 * np.arctan(lbda)
b = 0.5 * bscale * g
beta = np.sqrt(1.0 + b**2)
cos1 = np.cos(b * phi)
cos2 = np.cos(beta * phi)
sin1 = np.sin(b * phi)
sin2 = np.sin(beta * phi)
bg = np.sqrt(gamma_ref**2 - 1.0)

# exact transfer matrix
r11 = cos1 * (cos2 + lbda * sin2 / beta)
r12 = g * (1.0 + lbda**2) * cos1 * sin2 / beta
r13 = sin1 * (cos2 + lbda * sin2 / beta)
r14 = g * (1.0 + lbda**2) * sin1 * sin2 / beta
r21 = (
    cos1
    * (2.0 * beta * lbda * cos2 + (lbda**2 - beta**2) * sin2)
    / (g * beta * (1.0 + lbda**2))
)
r22 = cos1 * (cos2 + lbda * sin2 / beta)
r23 = (
    sin1
    * (2.0 * beta * lbda * cos2 + (lbda**2 - beta**2) * sin2)
    / (g * beta * (1.0 + lbda**2))
)
r24 = sin1 * (cos2 + lbda * sin2 / beta)
r31 = -r13
r32 = -r14
r33 = r22
r34 = r12
r41 = -r23
r42 = -r24
r43 = r21
r44 = r22
r55 = 1.0
r56 = 2.0 * lbda * g / bg**2
r66 = 1.0

# final data
xf = np.array([r11, r12, r13, r14, 0, 0])
pxf = np.array([r21, r22, r23, r24, 0, 0])
yf = np.array([r31, r32, r33, r34, 0, 0])
pyf = np.array([r41, r42, r43, r44, 0, 0])
tf = np.array([0, 0, 0, 0, 1.0, r56])
ptf = np.array([0, 0, 0, 0, 0, 1.0])

# compute differences
error_xi = np.max(np.abs(xi - initial["position_x"].to_numpy()))
error_pxi = np.max(np.abs(pxi - initial["momentum_x"].to_numpy()))
error_yi = np.max(np.abs(yi - initial["position_y"].to_numpy()))
error_pyi = np.max(np.abs(pyi - initial["momentum_y"].to_numpy()))
error_ti = np.max(np.abs(ti - initial["position_t"].to_numpy()))
error_pti = np.max(np.abs(pti - initial["momentum_t"].to_numpy()))

error_xf = np.max(np.abs(xf - final["position_x"].to_numpy()))
error_pxf = np.max(np.abs(pxf - final["momentum_x"].to_numpy()))
error_yf = np.max(np.abs(yf - final["position_y"].to_numpy()))
error_pyf = np.max(np.abs(pyf - final["momentum_y"].to_numpy()))
error_tf = np.max(np.abs(tf - final["position_t"].to_numpy()))
error_ptf = np.max(np.abs(ptf - final["momentum_t"].to_numpy()))

xf_max = np.max(np.abs(xf))
pxf_max = np.max(np.abs(pxf))
yf_max = np.max(np.abs(yf))
pyf_max = np.max(np.abs(pyf))
tf_max = np.max(np.abs(tf))
ptf_max = np.max(np.abs(ptf))

print("Initial maximum absolute difference in each column:")
print("Difference x, px, y, py, t, pt:")
print(error_xi)
print(error_pxi)
print(error_yi)
print(error_pyi)
print(error_ti)
print(error_pti)

atol = 1.0e-13
print(f"  atol={atol}")

assert np.allclose(
    [error_xi, error_pxi, error_yi, error_pyi, error_ti, error_pti],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    atol=atol,
)

print("")
print("Final maximum relative difference in each column:")
print("Difference x, px, y, py, t, pt:")
print(error_xf / xf_max)
print(error_pxf / pxf_max)
print(error_yf / yf_max)
print(error_pyf / pyf_max)
print(error_tf / tf_max)
print(error_ptf / ptf_max)

atol = 1.0e-7
print(f"  tol={atol}")

assert np.allclose(
    [
        error_xf / xf_max,
        error_pxf / pxf_max,
        error_yf / yf_max,
        error_pyf / pyf_max,
        error_tf / tf_max,
        error_ptf / ptf_max,
    ],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    atol=atol,
)
