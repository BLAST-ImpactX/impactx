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
alpha_x = rbc["alpha_x"]
beta_x = rbc["beta_x"]
alpha_y = rbc["alpha_y"]
beta_y = rbc["beta_y"]
dispersion_x = rbc["dispersion_x"]
dispersion_px = rbc["dispersion_px"]

sigma_xi = sigma_x.iloc[0]
sigma_yi = sigma_y.iloc[0]
sigma_ti = sigma_t.iloc[0]
emittance_xi = emittance_x.iloc[0]
emittance_yi = emittance_y.iloc[0]
emittance_ti = emittance_t.iloc[0]
alpha_xi = alpha_x.iloc[0]
beta_xi = beta_x.iloc[0]
alpha_yi = alpha_y.iloc[0]
beta_yi = beta_y.iloc[0]
dispersion_xi = dispersion_x.iloc[0]
dispersion_pxi = dispersion_px.iloc[0]

length = len(s) - 1

sf = s.iloc[length]
sigma_xf = sigma_x.iloc[length]
sigma_yf = sigma_y.iloc[length]
sigma_tf = sigma_t.iloc[length]
emittance_xf = emittance_x.iloc[length]
emittance_yf = emittance_y.iloc[length]
emittance_tf = emittance_t.iloc[length]
alpha_xf = alpha_x.iloc[length]
beta_xf = beta_x.iloc[length]
alpha_yf = alpha_y.iloc[length]
beta_yf = beta_y.iloc[length]
dispersion_xf = dispersion_x.iloc[length]
dispersion_pxf = dispersion_px.iloc[length]


print("Initial Beam:")
print(f"  sigx={sigma_xi:e} sigy={sigma_yi:e} sigt={sigma_ti:e}")
print(
    f"  emittance_x={emittance_xi:e} emittance_y={emittance_yi:e} emittance_t={emittance_ti:e}"
)

atol = 0.0  # ignored
rtol = 2.0e-2  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigma_xi, sigma_yi, sigma_ti, emittance_xi, emittance_yi, emittance_ti],
    [
        1.924711e-03,
        2.165646e-05,
        1.102534e-04,
        7.668809e-09,
        1.018986e-10,
        8.588054e-09,
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
rtol = 2.5e-2  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigma_xf, sigma_yf, sigma_tf, emittance_xf, emittance_yf, emittance_tf],
    [
        2.057123e-05,
        6.911405e-05,
        2.012573e-05,
        8.174618e-11,
        1.018986e-10,
        1.151058e-08,
    ],
    rtol=rtol,
    atol=atol,
)


print("")
print("Initial Twiss functions:")
print(
    f"  alpha_x={alpha_xi:e} beta_x={beta_xi:e} alpha_y={alpha_yi:e} beta_y={beta_yi:e}"
)
print(f"  dispersion_x={dispersion_xi:e} dispersion_px={dispersion_pxi:e}")

atol = 0.0  # ignored
rtol = 2.5e-2
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [alpha_xi, beta_xi, alpha_yi, beta_yi, dispersion_xi],
    [
        1.340770e00,
        1.440253e01,
        -1.347747e00,
        4.602637e00,
        -2.667038e-01,
    ],
    rtol=rtol,
    atol=atol,
)

print("")
print("Final Twiss functions:")
print(
    f"  alpha_x={alpha_xf:e} beta_x={beta_xf:e} alpha_y={alpha_yf:e} beta_y={beta_yf:e}"
)
print(f"  dispersion_x={dispersion_xf:e} dispersion_px={dispersion_pxf:e}")

atol = 0.0  # ignored
rtol = 2.0e-2
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [beta_xf, alpha_yf, beta_yf],
    [
        5.176642e00,
        -4.973010e00,
        4.687750e01,
    ],
    rtol=rtol,
    atol=atol,
)

# We use absolute tolerance for the following quantities, because these

assert np.allclose(
    [alpha_xf],
    [
        0.0,
    ],
    atol=0.1,
)
assert np.allclose(
    [dispersion_xf],
    [
        0.0,
    ],
    atol=4.0e-5,
)
