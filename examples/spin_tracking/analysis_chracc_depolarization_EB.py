#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import glob
import re
from pathlib import Path

import numpy as np
import openpmd_api as io
import pandas as pd
from scipy import special


# Load data from envelope simulation
def read_time_series(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.

    Returns
    -------
    pandas.DataFrame
    """

    def read_file(file_pattern):
        for filename in glob.glob(file_pattern):
            df = pd.read_csv(filename, delimiter=r"\s+")
            if "step" not in df.columns:
                step = int(re.findall(r"[0-9]+", filename)[0])
                df["step"] = step
            yield df

    return pd.concat(
        read_file(file_pattern),
        axis=0,
        ignore_index=True,
    )  # .set_index('id')


def data_is_double(file_pattern):
    """Detect float precision (single vs double) from a text diagnostic's digits."""
    text = "".join(Path(f).read_text() for f in glob.glob(file_pattern))
    return any(
        len(m.replace(".", "") if e else m.replace(".", "").strip("0")) >= 12
        for m, e in re.findall(r"(\d*\.\d+)([eE][+-]?\d+)?", text)
    )


# collect reduced beam characteristics
rbc = read_time_series("diags/reduced_beam_characteristics.*")
is_double = data_is_double("diags/reduced_beam_characteristics.*")

# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]

beam_initial = series.iterations[1].particles["beam"]
gamma_i = beam_initial.get_attribute("gamma_ref")
beta_gamma_i = beam_initial.get_attribute("beta_gamma_ref")
beam_final = series.iterations[last_step].particles["beam"]
gamma_f = beam_final.get_attribute("gamma_ref")
beta_gamma_f = beam_final.get_attribute("beta_gamma_ref")

# numerical parameters based on input file
gyro_anomaly = 1.7928473446  # for protons
ez_value = 10.0  # accelerating field in 1/m (dgamma/dz) from input
bz_value = 3.0  # longitudinal magnetic field in 1/m from input
ds_value = 50.0  # accelerating distance in m from input
sigmaX = 0.003  # value lambdaX from input
sigmaPx = 0.2  # value lambdaPx from input

Pxi = 0.4  # polarization_x from input
Pyi = 0.8  # polarization_y from input
Pzi = 0.2  # polarization_z from input

print("Initial Beam:")
polarization_x = rbc["mean_sx"].iloc[0]
polarization_y = rbc["mean_sy"].iloc[0]
polarization_z = rbc["mean_sz"].iloc[0]
print(
    f"  polarization_x={polarization_x:e} polarization_y={polarization_y:e} polarization_z={polarization_z:e}"
)

# Number of particles (from input file)
num_particles = 100000

atol = 1.3 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  atol={atol}")

assert np.allclose(
    [polarization_x, polarization_y, polarization_z],
    [
        Pxi,
        Pyi,
        Pzi,
    ],
    rtol=0.0,
    atol=atol,
)

# predicted final polarization
alpha = bz_value / 2.0
sigmaPx_dyn = sigmaPx * beta_gamma_i
size_factor = np.sqrt(sigmaPx_dyn**2 + alpha**2 * sigmaX**2)
theta = alpha / ez_value * np.log((gamma_f + beta_gamma_f) / (gamma_i + beta_gamma_i))
mu = 0.5 * (beta_gamma_f / (1 + gamma_f) - beta_gamma_i / (1 + gamma_i))
sinG = np.sin(2.0 * gyro_anomaly * theta)
cosG = np.cos(2.0 * gyro_anomaly * theta)
F_func = (
    beta_gamma_i / (1 + gamma_i) * sinG
    + ez_value * (1 - cosG) / (2.0 * alpha)
    + 2 * mu * sinG
)
G_func = (
    beta_gamma_i / (1 + gamma_i) * (1 - cosG)
    - ez_value * sinG / (2.0 * alpha)
    - 2 * mu * cosG
)
damping_eigenvalue = np.sqrt(F_func**2 + G_func**2) * size_factor

# DawsonF damping factor
ev_arg = damping_eigenvalue / np.sqrt(2.0)
dawson_factor = ev_arg * special.dawsn(ev_arg)

# Factor due to mixing
damping_factor_xy = 1 - dawson_factor
damping_factor_z = 1 - 2.0 * dawson_factor
Px_damped = damping_factor_xy * Pxi
Py_damped = damping_factor_xy * Pyi
Pz_damped = damping_factor_z * Pzi

# Design rotation
ref_spin_angle = 2.0 * theta * (1 + gyro_anomaly)
Pxf = Px_damped * np.cos(ref_spin_angle) + Py_damped * np.sin(ref_spin_angle)
Pyf = -Px_damped * np.sin(ref_spin_angle) + Py_damped * np.cos(ref_spin_angle)
Pzf = Pz_damped

print("")
print("Final Predicted Polarization:")
print(f"  polarization_x={Pxf:e} polarization_y={Pyf:e} polarization_z={Pzf:e}")

print("")
print("Final Beam:")
polarization_x = rbc["mean_sx"].iloc[-1]
polarization_y = rbc["mean_sy"].iloc[-1]
polarization_z = rbc["mean_sz"].iloc[-1]
print(
    f"  polarization_x={polarization_x:e} polarization_y={polarization_y:e} polarization_z={polarization_z:e}"
)

atol = 1.5 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  atol={atol}")

assert np.allclose(
    [polarization_x, polarization_y, polarization_z],
    [
        Pxf,
        Pyf,
        Pzf,
    ],
    rtol=0.0,
    atol=atol,
)

# numerical tests of spin moment conditions
sigma_sx = rbc["sigma_sx"].iloc[-1]
sigma_sy = rbc["sigma_sy"].iloc[-1]
sigma_sz = rbc["sigma_sz"].iloc[-1]
polarization = np.sqrt(polarization_x**2 + polarization_y**2 + polarization_z**2)
condition = sigma_sx**2 + sigma_sy**2 + sigma_sz**2 + polarization**2

print("")
print(f"Spin moment consistency condition = {condition:e}")

atol = 1.0e-12 if is_double else 1.0e-3
print(f"  atol={atol}")

assert np.allclose(
    [condition],
    [
        1.0,
    ],
    rtol=0.0,
    atol=atol,
)
