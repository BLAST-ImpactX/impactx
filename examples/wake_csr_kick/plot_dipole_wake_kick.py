#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse

import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io
import scipy.constants as sc

# options to run this script
parser = argparse.ArgumentParser(
    description="Plot the wake-induced energy kick for a bunch propagating in an RF cavity."
)
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

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
mp_eV = mp * clite**2 / qe

# Basic beam parameters
bunch_charge_C = 1.3252572e-8
kin_energy_MeV = 800.0
gamma = 1.0e6 * kin_energy_MeV / mp_eV + 1.0
beta_gamma_2 = gamma**2 - 1.0
beta = np.sqrt(beta_gamma_2) / gamma
p0_eVc = mp_eV * beta * gamma

# Beam distribution parameters
sigmact = 1.0e-3
sigmaz = beta * sigmact

# Simulation slice length
L = 2.889612

# Short-range wake intensity parameter
K0 = 4.30355e13  # first fitted parameter
s0 = 0.00175  # second fitted parameter (in meters)
Kscale_eV_per_m = K0 * bunch_charge_C

# Initial particle data

ti = initial_sort["position_t"]
zi = -beta * ti
zi_scaled = zi / sigmaz

# Predicted wake for a uniform longitudinal distribution (KVdist):
a = sigmaz * np.sqrt(3.0)
radical = np.sqrt((a - zi) / s0)
exp_plus = np.exp(radical)
exp_minus = np.exp(-radical)
wake = -Kscale_eV_per_m * exp_minus * (s0 - exp_plus * s0 + radical * s0) / a

# Predicted particle kick: the wake decelerates
dE_predicted = -wake * L / p0_eVc
pt_predicted = -dE_predicted

# Maximum momentum kick
pt_max = pt_predicted.abs().max()

# Final particle data

ptf = final_sort["momentum_t"]
dEf = -ptf

# Difference between value and prediction


dpt_max = (ptf - pt_predicted).abs().max()
dpt_rms = np.sqrt(np.mean(np.square(ptf - pt_predicted)))

print()
print("Difference between predicted and computed final momentum, absolute rms:")
print("dpt_rms", dpt_rms)

print()
print("Difference between predicted and computed final momentum, absolute max:")
print("dpt_max", dpt_max)

print()
print("Difference between predicted and computed final momentum (rms), relative:")
print("dpt_rms/pt_max", dpt_rms / pt_max)

print()
print("Difference between predicted and computed final momentum (max), relative:")
print("dpt_max/pt_max", dpt_max / pt_max)

# Plotting
plt.figure(figsize=(10, 6))
plt.xlabel(r"initial $z/\sigma_z$", fontsize=25)
plt.ylabel(r"$\Delta E/cp_0$", fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(True)

# Plot the data
plt.plot(zi_scaled, dEf, "bo", label="computed")
plt.plot(zi_scaled, dE_predicted, "ro", markersize=5, label="predicted")
plt.legend(fontsize=20)
plt.tight_layout()
if args.save_png:
    plt.savefig("User_wake_longitudinal_kick.png")
else:
    plt.show()
