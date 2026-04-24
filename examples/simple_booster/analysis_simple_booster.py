#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import glob
import math
import re

import openpmd_api as io
import pandas as pd
from booster_impactx_lattice import get_lattice

from impactx import elements


def read_file(file_pattern):
    for filename in glob.glob(file_pattern):
        df = pd.read_csv(filename, delimiter=r"\s+")
        if "step" not in df.columns:
            step = int(re.findall(r"[0-9]+", filename)[0])
            df["step"] = step
        else:
            df = df[df["step"] != "step"]
        df = df.apply(pd.to_numeric)
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


# load the Booster lattice
lattice = elements.KnownElementsList()
lattice.extend(get_lattice())

# total lattice length (nominal Fermilab Booster circumference: 474.20 m)
total_s = sum(element.ds for element in lattice)
expected_s = 474.20
# tolerance loose enough to hold in single precision accumulation across the ring
assert math.isclose(total_s, expected_s, rel_tol=1.0e-4), (
    f"lattice length {total_s} m does not match expected {expected_s} m"
)

series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

beta_ref = series.iterations[1].particles["beam"].get_attribute("beta_ref")

rbc = read_time_series("diags/reduced_beam_characteristics.*")

sigma_x = rbc["sigma_x"]
sigma_y = rbc["sigma_y"]
sigma_t = rbc["sigma_t"]

sigma_px = rbc["sigma_px"]
sigma_py = rbc["sigma_py"]
sigma_pt = rbc["sigma_pt"]

emittance_x = rbc["emittance_x"]
emittance_y = rbc["emittance_y"]
emittance_t = rbc["emittance_t"]

print("\t\tInitial\t\tFinal")
print("\t\t=======\t\t=====")
print(f"std x\t\t{sigma_x.iloc[0]:7g}\t{sigma_x.iloc[-1]:7g}")
print(f"std px\t\t{sigma_px.iloc[0]:7g}\t{sigma_px.iloc[-1]:7g}")
print(f"emittance x\t{emittance_x.iloc[0]:7g}\t{emittance_x.iloc[-1]:7g}")
print()
print(f"std y\t\t{sigma_y.iloc[0]:7g}\t{sigma_y.iloc[-1]:7g}")
print(f"std py\t\t{sigma_py.iloc[0]:7g}\t{sigma_py.iloc[-1]:7g}")
print(f"emittance y\t{emittance_y.iloc[0]:7g}\t{emittance_y.iloc[-1]:7g}")
print()
print(f"std t\t\t{sigma_t.iloc[0]:7g}\t\t{sigma_t.iloc[-1]:7g}")
print(f"std pt\t\t{sigma_pt.iloc[0]:7g}\t{sigma_pt.iloc[-1]:7g}")
print(f"emittance t\t{emittance_t.iloc[0]:7g}\t{emittance_t.iloc[-1]:7g}")

assert math.isclose(sigma_x.iloc[0], 0.0083783, rel_tol=1.0e-4)
assert math.isclose(sigma_x.iloc[-1], 0.0083783, rel_tol=1.0e-2), (
    "sigma_x change over one turn too large"
)

assert math.isclose(sigma_px.iloc[0], 0.000226298, rel_tol=1.0e-4)
assert math.isclose(sigma_px.iloc[-1], 0.000226298, rel_tol=1.0e-2), (
    "sigma_px change over one turn too large"
)

assert math.isclose(emittance_x.iloc[0], 1.89474e-06, rel_tol=1.0e-4)
assert math.isclose(emittance_x.iloc[-1], 1.89474e-06, rel_tol=1.0e-2), (
    "emittance_x change over one turn too large"
)

assert math.isclose(sigma_y.iloc[0], 0.00299246, rel_tol=1.0e-4)
assert math.isclose(sigma_y.iloc[-1], 0.00299246, rel_tol=1.0e-2), (
    "sigma_y change over one turn too large"
)

assert math.isclose(sigma_py.iloc[0], 0.000575451, rel_tol=1.0e-4)
assert math.isclose(sigma_py.iloc[-1], 0.000575451, rel_tol=1.0e-2), (
    "sigma_py change over one turn too large"
)

assert math.isclose(emittance_y.iloc[0], 1.72198e-06, rel_tol=1.0e-4)
assert math.isclose(emittance_y.iloc[-1], 1.72198e-06, rel_tol=1.0e-2), (
    "emittance_y change over one turn too large"
)

assert math.isclose(sigma_t.iloc[0], 1.16878, rel_tol=1.0e-4)
assert math.isclose(sigma_t.iloc[-1], 1.16878, rel_tol=1.0e-2), (
    "sigma_t change over one turn too large"
)

assert math.isclose(sigma_pt.iloc[0], 0.000914371, rel_tol=1.0e-4)
assert math.isclose(sigma_pt.iloc[-1], 0.000914371, rel_tol=1.0e-2), (
    "sigma_pt change over one turn too large"
)

assert math.isclose(emittance_t.iloc[0], 0.00106837, rel_tol=1.0e-4)
assert math.isclose(emittance_t.iloc[-1], 0.00106837, rel_tol=1.0e-2), (
    "emittance_t change over one turn too large"
)
