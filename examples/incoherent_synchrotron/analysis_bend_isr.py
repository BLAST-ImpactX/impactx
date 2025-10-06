#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import glob
import re

import numpy as np
import pandas as pd


def read_file(file_pattern):
    for filename in glob.glob(file_pattern):
        df = pd.read_csv(filename, delimiter=r"\s+")
        if "step" not in df.columns:
            step = int(re.findall(r"[0-9]+", filename)[0])
            df["step"] = step
        yield df


def read_time_series(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.

    Returns
    -------
    pandas.DataFrame
    """
    return pd.concat(
        read_file(file_pattern),
        axis=0,
        ignore_index=True,
    )  # .set_index('id')


# read reduced diagnostics
rbc = read_time_series("diags/reduced_beam_characteristics.*")

s = rbc["s"]
px_mean = rbc["mean_px"]
py_mean = rbc["mean_py"]
pt_mean = rbc["mean_pt"]
sigma_px = rbc["sigma_px"]
sigma_py = rbc["sigma_py"]
sigma_pt = rbc["sigma_pt"]

px_meani = px_mean.iloc[0]
py_meani = py_mean.iloc[0]
pt_meani = pt_mean.iloc[0]
sig_pxi = sigma_px.iloc[0]
sig_pyi = sigma_py.iloc[0]
sig_pti = sigma_pt.iloc[0]

length = len(s) - 1

sf = s.iloc[length]

px_meanf = px_mean.iloc[length]
py_meanf = py_mean.iloc[length]
pt_meanf = pt_mean.iloc[length]
sig_pxf = sigma_px.iloc[length]
sig_pyf = sigma_py.iloc[length]
sig_ptf = sigma_pt.iloc[length]

print("Initial Beam:")
print(f"  px_mean={px_meani:e} py_mean={py_meani:e} pt_mean={pt_meani:e}")
print(f"  sig_px={sig_pxi:e} sig_py={sig_pyi:e} sig_pt={sig_pti:e}")


atol = 1.0e-6
print(f"  atol={atol}")
assert np.allclose(
    [px_meani, py_meani, pt_meani, sig_pxi, sig_pyi, sig_pti],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    atol=atol,
)

# Physical constants:
re_classical = 2.8179403205e-15  # classical electron radius
lambda_compton_reduced = 3.8615926744e-13  # reduced Compton wavelength

# Problem parameters:
ds = 0.1
rc = 10.0
gamma = 195696.117901
num_particles = 10000

# Characteristic length of energy loss:
length_isr = 3.0 / 2.0 * rc**2 / (re_classical * gamma**3)

print("")
print("Length and characteristic length for energy loss [m]:")
print(f" Length={ds:e} Length_ISR={length_isr:e}")

# Predicted energy loss and energy spread:
dpt = 2.0 / 3.0 * re_classical * ds / (rc**2) * gamma**3

dsigpt2 = (
    55.0
    / (24 * np.sqrt(3.0))
    * lambda_compton_reduced
    * re_classical
    * ds
    / rc**3
    * gamma**5
)
dsigpt = np.sqrt(dsigpt2)

print("")
print("Final Beam:")
print(f" pt_mean={pt_meanf:e} sig_pt={sig_ptf}")
print("Predicted (assuming that Length << Length_isr):")
print(f" pt_mean={dpt:e} sig_pt={dsigpt:e}")
print("")

rtol = 10.0 * num_particles**-0.5  # from random sampling of a smooth distribution
assert np.allclose(
    [pt_meanf, sig_ptf],
    [
        dpt,
        dsigpt,
    ],
    rtol=rtol,
    atol=atol,
)
