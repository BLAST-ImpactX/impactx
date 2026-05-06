#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import argparse
import glob
import re

import matplotlib.pyplot as plt
import openpmd_api as io
import pandas as pd
import scipy


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


# options to run this script
parser = argparse.ArgumentParser(description="Plot the FODO benchmark.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
first_step = list(series.iterations)[0]
initial = series.iterations[first_step].particles["beam"].to_df()
last_step = list(series.iterations)[-1]
final = series.iterations[last_step].particles["beam"].to_df()

beta_ref = series.iterations[first_step].particles["beam"].get_attribute("beta_ref")

# columns in rbc file
# step s mean_x min_x max_x mean_y min_y max_y mean_t min_t max_t sigma_x sigma_y sigma_t mean_px min_px max_px mean_py min_py max_py mean_pt min_pt max_pt sigma_px sigma_py sigma_pt emittance_x emittance_y emittance_t alpha_x alpha_y alpha_t beta_x beta_y beta_t dispersion_x dispersion_px dispersion_y dispersion_py emittance_xn emittance_yn emittance_tn charge_C

rbc = read_time_series("diags/reduced_beam_characteristics.*")

s = rbc["s"]

min_x = rbc["min_x"]
max_x = rbc["max_x"]

min_y = rbc["min_y"]
max_y = rbc["max_y"]

min_t = rbc["min_t"]
max_t = rbc["max_t"]

sigma_x = rbc["sigma_x"]
sigma_y = rbc["sigma_y"]
sigma_t = rbc["sigma_t"]
sigma_pt = rbc["sigma_pt"]

charge = rbc["charge_C"] / scipy.constants.eV

f, ax = plt.subplots(2, 2, sharex="row")

ax[0, 0].set_title(r"$\sigma_x$")
ax[0, 0].plot(s, sigma_x * 1000.0)

ax[0, 0].set_ylabel(r"$\sigma_x$ [mm]")

ax[0, 1].set_title(r"$\sigma_y$")
ax[0, 1].plot(s, sigma_y * 1000.0)
ax[0, 1].set_ylabel(r"$\sigma_y$ [mm]")

ax[1, 0].set_title(r"$\sigma_t$")
ax[1, 0].plot(s, sigma_t, label=r"$sigma_t$ [m]")
ax[1, 0].set_xlabel("s [m]")
ax[1, 0].set_ylabel(r"$\sigma_t$ [m]")

ax[1, 1].set_title(r"$\sigma_{pt}$")
ax[1, 1].plot(s, sigma_pt * 1000.0, label="sigma_pt")
ax[1, 1].set_ylabel(r"$\sigma_{pt} \, [\times 1000]$")
ax[1, 1].set_xlabel("s [m]")

plt.tight_layout()

if args.save_png:
    plt.savefig("simple_booster_sigma.png")
else:
    plt.show()

f, ax = plt.subplots(2, 2, sharex="row", sharey="row")

ax[0, 0].set_title("x vs. px initial")
ax[0, 0].plot(
    initial["position_x"] * 1000, initial["momentum_x"] * 1000, ".", label="initial x vs. px",
)
ax[0, 0].set_xlabel("x [mm]")
ax[0, 0].set_ylabel("px [x1000]")

ax[0, 1].set_title("x vs. px final")
ax[0, 1].plot(
    final["position_x"] * 1000, final["momentum_x"] * 1000, ".", label="final x vs. px",
)
ax[0, 1].set_xlabel("x [mm]")

ax[1, 0].set_title("y vs. py initial")
ax[1, 0].plot(
    initial["position_y"] * 1000, initial["momentum_y"] * 1000, ".", label="initial y vs. py",
)
ax[1, 0].set_xlabel("y [mm]")
ax[1, 0].set_ylabel("py [x1000]")
# ax[1, 0].legend(loc="best")

ax[1, 1].set_title("y vs. py final")
ax[1, 1].plot(
    final["position_y"] * 1000, final["momentum_y"] * 1000, ".", label="final y vs. py",
)
ax[1, 1].set_xlabel("y [mm]")


plt.tight_layout()

if args.save_png:
    plt.savefig("simple_booster_scatter.png")
else:
    plt.show()
