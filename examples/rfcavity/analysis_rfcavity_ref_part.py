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


# read reference particle data
rbc = read_time_series("diags/ref_particle.*")

s = rbc["s"]
gamma = rbc["gamma"]

si = s.iloc[0]
gammai = gamma.iloc[0]

sf = s.iloc[-1]
gammaf = gamma.iloc[-1]

print("")
print("Initial Beam:")
print(f"  s_ref={si:e} gamma_ref={gammai:e}")

atol = 1.0e-4  # ignored
print(f"  atol={atol}")

assert np.allclose(
    [si, gammai],
    [
        0.000000,
        451.09877160930125,
    ],
    atol=atol,
)


print("")
print("Final Beam:")
print(f"  s_ref={sf:e} gamma_ref={gammaf:e}")

atol = 1.0e-4  # ignored
print(f"  atol={atol}")

assert np.allclose(
    [sf, gammaf],
    [
        5.9391682799999987,
        585.11989430128892,
    ],
    atol=atol,
)
