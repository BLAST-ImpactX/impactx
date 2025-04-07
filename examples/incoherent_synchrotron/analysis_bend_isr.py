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
px_mean = rbc["px_mean"]
py_mean = rbc["py_mean"]
pt_mean = rbc["pt_mean"]
sig_px = rbc["sig_px"]
sig_py = rbc["sig_py"]
sig_pt = rbc["sig_pt"]

px_meani = px_mean.iloc[0]
py_meani = py_mean.iloc[0]
pt_meani = pt_mean.iloc[0]
sig_pxi = sig_px.iloc[0]
sig_pyi = sig_py.iloc[0]
sig_pti = sig_pt.iloc[0]

length = len(s) - 1

sf = s.iloc[length]

px_meanf = px_mean.iloc[length]
py_meanf = py_mean.iloc[length]
pt_meanf = pt_mean.iloc[length]
sig_pxf = sig_px.iloc[length]
sig_pyf = sig_py.iloc[length]
sig_ptf = sig_pt.iloc[length]

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
num_particles = 100000

# Predicted energy loss and energy spread:
dpt = 2.0 / 3.0 * re_classical * ds / (rc**2) * gamma**3

dsigpt2 = (
    55.0 / (24 * np.sqrt(3.0)) * lambda_compton_reduced * re_classical * ds / rc**3 * gamma**5
)
dsigpt = np.sqrt(dsigpt2)

print("")
print("Final Beam:")
print(f" pt_mean={pt_meanf:e} sig_pt={sig_ptf}")
print("Predicted:")
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
