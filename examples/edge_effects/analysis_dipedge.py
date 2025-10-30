#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import math

import numpy as np
import openpmd_api as io

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# Basic input parameters
g = 1.0e-3
phi = math.pi / 8.0
rc = 10.0
R = 1.0
K0 = math.pi**2 / 6.0
K3 = 1.0 / 6.0
Kar = [K0, 0, 0, K3, 0, 0, 0]
delta = 0.0

# Derived quantities
cs = math.cos(phi)
sn = math.sin(phi)
tn = sn / cs
sc = 1.0 / cs

# Lie generator coefficients
c1 = g * Kar[1] / (rc * cs)
c2 = sn * g**2 * Kar[0] / rc**2 * 1.0 / (2.0 * cs**3 * (1 + delta))
c3 = g**2 / rc * Kar[0] / (cs**2 * (1 + delta))
c4 = 1 / (1 + delta) * g / rc * Kar[1] * sn / cs**2
c5 = sn / cs * 1.0 / (2 * rc)
c6 = g * Kar[1] / rc * sn**2 / (4 * rc * (1 + delta) * cs**3)
c7 = (
    1
    / (2 * cs**3 * (1 + delta))
    * (g * Kar[1] / (2 * rc**2) + (1 + sn**2) * g / rc**2 * Kar[2])
)
c8 = 1 / 6 * tn**3 / (2 * rc**2 * (1 + delta))
c9 = 1 / 2 * (tn * sc**2 / (2 * rc**2 * (1 + delta)))
c10 = 1 / (2 * (1 + delta)) * tn**2 / rc
c11 = 1 / (2 * rc * (1 + delta))
c12 = 1 / 24 * (4 / cs - 8 / cs**3) * Kar[3] / (rc**2 * g * (1 + delta))
c13 = sn**2 / (2 * cs**3) * g**2 / (rc * R) * Kar[4]
c14 = 1 / 2 * sn / cs**3 * g / (rc * R) * Kar[5]
c15 = Kar[6] / (rc * R) * 1 / cs**3

xi = initial["position_x"]
pxi = initial["momentum_x"]
yi = initial["position_y"]
pyi = initial["momentum_y"]
ti = initial["position_t"]
pti = initial["momentum_t"]

Omega_initial = (
    xi * c1
    - xi * c2
    + pxi * c3
    + (xi * pxi - yi * pyi) * c4
    + (xi**2 - yi**2) * c5
    - xi**2 * c6
    + yi**2 * c7
    - xi**3 * c8
    + xi * yi**2 * c9
    + (xi**2 * pxi - yi**2 * pxi - 2 * xi * yi * pyi) * c10
    - yi**2 * pxi * c11
    + yi**4 * c12
    + xi * c13
    + (yi**2 - xi**2) * c14
    + (xi * yi**2 / 2 - xi**3 / 6) * c15
)

xf = final["position_x"]
pxf = final["momentum_x"]
yf = final["position_y"]
pyf = final["momentum_y"]
tf = final["position_t"]
ptf = final["momentum_t"]

Omega_final = (
    xf * c1
    - xf * c2
    + pxf * c3
    + (xf * pxf - yf * pyf) * c4
    + (xf**2 - yf**2) * c5
    - xf**2 * c6
    + yf**2 * c7
    - xf**3 * c8
    + xf * yf**2 * c9
    + (xf**2 * pxf - yf**2 * pxf - 2 * xf * yf * pyf) * c10
    - yf**2 * pxf * c11
    + yf**4 * c12
    + xf * c13
    + (yf**2 - xf**2) * c14
    + (xf * yf**2 / 2 - xf**3 / 6) * c15
)

Delta_Omega = (Omega_final - Omega_initial).abs()

dx = (xf - xi).abs().max()
dpx = (pxf - pxi).abs().max()
dy = (yf - yi).abs().max()
dpy = (pyf - pyi).abs().max()
dt = (tf - ti).abs().max()
dpt = (ptf - pti).abs().max()

print("Change in the coordinates and momenta:")
print("dx", dx)
print("dpx", dpx)
print("dy", dy)
print("dpy", dpy)
print("dt", dt)
print("dpt", dpt)

print("Change in the Lie generator, for each initial condition:")
print(Delta_Omega)

atol = 1.5e-10
print(f"  atol={atol}")

assert np.allclose(
    [Delta_Omega.max()],
    [
        0.0,
    ],
    atol=atol,
)
