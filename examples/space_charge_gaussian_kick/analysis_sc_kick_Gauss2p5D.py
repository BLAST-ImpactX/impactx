#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
import scipy.constants as sc

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# Physical constants
qe = sc.elementary_charge
me = sc.electron_mass
clite = sc.speed_of_light
eps0 = sc.epsilon_0
me_eV = me * clite**2 / qe

# Basic beam parameters
bunch_charge_C = 1.0e-9
kin_energy_MeV = 100
gamma = 1.0e6 * kin_energy_MeV / me_eV + 1.0
beta_gamma_2 = gamma**2 - 1.0
beta = np.sqrt(beta_gamma_2) / gamma

# Space charge intensity parameter (units of length, meters)
Kscale = (
    qe * bunch_charge_C / (4.0 * np.pi * eps0 * me * clite**2 * beta_gamma_2 * gamma)
)

# Beam distribution parameters
sigmax = 4.0e-5
sigmay = 4.0e-5
sigmact = 1.0e-3
sigmaz = beta * sigmact

# Simulation slice length
L = 1.0

# Initial particle data

xi = initial["position_x"]
pxi = initial["momentum_x"]
yi = initial["position_y"]
pyi = initial["momentum_y"]
ti = initial["position_t"]
pti = initial["momentum_t"]
zi = beta * ti

# Predicted momentum kick, from J. Qiang, Phys. Rev. Accel. Beams 28, 114602 (2025), eqs. (31-32)

ri_2 = xi**2 + yi**2
gauss_exp = np.exp(-ri_2 / (2.0 * sigmax**2))
gauss_xy_factor = 2.0 * (1.0 - gauss_exp) / ri_2
gauss_pdf_z = np.exp(-(zi**2) / (2.0 * sigmaz**2)) / (np.sqrt(2.0 * np.pi) * sigmaz)
px_predicted = L * Kscale * gauss_pdf_z * xi * gauss_xy_factor
py_predicted = L * Kscale * gauss_pdf_z * yi * gauss_xy_factor

# Maximum momentum kick
px_max = px_predicted.abs().max()
py_max = py_predicted.abs().max()

# Final particle data

xf = final["position_x"]
pxf = final["momentum_x"]
yf = final["position_y"]
pyf = final["momentum_y"]
tf = final["position_t"]
ptf = final["momentum_t"]

# Difference between value and prediction

dpx = (pxf - px_predicted).abs().max()
dpy = (pyf - py_predicted).abs().max()
# dpt = (ptf - pti).abs().max()

print("Difference between predicted and computed final momentum, absolute:")
print("dpx", dpx)
print("dpy", dpy)
# print("dpt", dpt)

print("Difference between predicted and computed final momentum, relative:")
print("dpx/px_max", dpx / px_max)
print("dpy/py_max", dpy / py_max)
# print("dpt", dpt)

atol = 5.1e-2
print(f"  tol={atol}")

assert np.allclose(
    [dpx / px_max, dpy / py_max],
    [0.0, 0.0],
    atol=atol,
)
