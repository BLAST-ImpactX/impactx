#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
import scipy.constants as sc
import pandas as pd

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
beam_initial = series.iterations[1].particles["beam"]
initial_sort = beam_initial.to_df().set_index("id")
beam_final = series.iterations[last_step].particles["beam"]
final_sort = beam_final.to_df().set_index("id")

# Physical constants
qe = sc.elementary_charge
mp = sc.proton_mass
clite = sc.speed_of_light
eps0 = sc.epsilon_0
mp_eV = mp * clite**2 / qe

# Basic beam parameters
R0 = 1.0e-3  # initial beam radius (m)
Ib = 0.15  # beam current (A)
kin_energy_MeV = 250  # reference energy
gamma = 1.0e6 * kin_energy_MeV / mp_eV + 1.0
beta_gamma_3 = (np.sqrt(gamma**2 - 1.0)) ** 3

# Alfven current
IA = (4.0 * np.pi * eps0 * mp * clite**3) / qe
print("IA = ")
print(IA)

# Generalized perveance
Kpv = 2.0 * Ib / (beta_gamma_3 * IA)
print("Kpv = ")
print(Kpv)

# Simulation slice length
double_constant = 1.516770632602484  # constant independent of beam parameters
doubling_distance = R0 * double_constant / np.sqrt(Kpv)
print("doubling distance = ")
print(doubling_distance)

# Initial particle data

xi = initial_sort["position_x"]
pxi = initial_sort["momentum_x"]
yi = initial_sort["position_y"]
pyi = initial_sort["momentum_y"]

# Predicted momentum kick

px_predicted = doubling_distance * Kpv * xi / R0**2
py_predicted = doubling_distance * Kpv * yi / R0**2

# Maximum predicted momentum kick
pr_max = doubling_distance * Kpv / R0

print()
print("Maximum predicted momentum kick:")
print("pr_max", pr_max)

# Final particle data

xf = final_sort["position_x"]
pxf = final_sort["momentum_x"]
yf = final_sort["position_y"]
pyf = final_sort["momentum_y"]

# Difference between value and prediction

dpx_max = (pxf - px_predicted).abs().max()
dpy_max = (pyf - py_predicted).abs().max()

dpx_rms = np.sqrt(np.mean(np.square(pxf - px_predicted)))
dpy_rms = np.sqrt(np.mean(np.square(pyf - py_predicted)))

print()
print("Difference between predicted and computed final momentum, absolute rms:")
print("dpx_rms", dpx_rms)
print("dpy_rms", dpy_rms)

print()
print("Difference between predicted and computed final momentum, absolute max:")
print("dpx_max", dpx_max)
print("dpy_max", dpy_max)

print()
print("Difference between predicted and computed final momentum (rms), relative:")
print("dpx_rms/pr_max", dpx_rms / pr_max)
print("dpy_rms/pr_max", dpy_rms / pr_max)

# Test maximum error:
atol = 5.1e-2
print(f"  tol={atol}")

assert np.allclose(
    [dpx_rms / pr_max, dpy_rms / pr_max],
    [0.0, 0.0],
    atol=atol,
)

print()
print("Difference between predicted and computed final momentum (max), relative:")
print("dpx_max/pr_max", dpx_max / pr_max)
print("dpy_max/pr_max", dpy_max / pr_max)

# Test maximum error:
atol = 6.1e-2
print(f"  tol={atol}")

assert np.allclose(
    [dpx_max / pr_max, dpy_max / pr_max],
    [0.0, 0.0],
    atol=atol,
)
