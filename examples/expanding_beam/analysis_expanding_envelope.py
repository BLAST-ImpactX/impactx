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
sigma_x = rbc["sigma_x"]
sigma_y = rbc["sigma_y"]
sigma_t = rbc["sigma_t"]
emittance_x = rbc["emittance_x"]
emittance_y = rbc["emittance_y"]
emittance_t = rbc["emittance_t"]

sigma_xi = sigma_x.iloc[0]
sigma_yi = sigma_y.iloc[0]
sigma_ti = sigma_t.iloc[0]
emittance_xi = emittance_x.iloc[0]
emittance_yi = emittance_y.iloc[0]
emittance_ti = emittance_t.iloc[0]

length = len(s) - 1

sf = s.iloc[length]
sigma_xf = sigma_x.iloc[length]
sigma_yf = sigma_y.iloc[length]
sigma_tf = sigma_t.iloc[length]
emittance_xf = emittance_x.iloc[length]
emittance_yf = emittance_y.iloc[length]
emittance_tf = emittance_t.iloc[length]


print("Initial Beam:")
print(f"  sigx={sigma_xi:e} sigy={sigma_yi:e} sigt={sigma_ti:e}")
print(
    f"  emittance_x={emittance_xi:e} emittance_y={emittance_yi:e} emittance_t={emittance_ti:e}"
)

atol = 0.0  # ignored
rtol = 1.0e-3
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigma_xi, sigma_yi, sigma_ti, emittance_xi, emittance_yi, emittance_ti],
    [
        4.4721359550e-004,
        4.4721359550e-004,
        9.1224186858e-007,
        0.0e-006,
        0.0e-006,
        0.0e-006,
    ],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam:")
print(f"  sigx={sigma_xf:e} sigy={sigma_yf:e} sigt={sigma_tf:e}")
print(
    f"  emittance_x={emittance_xf:e} emittance_y={emittance_yf:e} emittance_t={emittance_tf:e}"
)

atol = 0.0  # ignored
rtol = 1.0e-3
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigma_xf, sigma_yf, sigma_tf],
    [
        9.029112e-04,
        9.029112e-04,
        1.841402e-06,
    ],
    rtol=rtol,
    atol=atol,
)
