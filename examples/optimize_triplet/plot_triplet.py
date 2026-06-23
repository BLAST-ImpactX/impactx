#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse
import glob
import re

import matplotlib.pyplot as plt
import pandas as pd

# options to run this script
parser = argparse.ArgumentParser(description="Plot the quadrupole triplet benchmark.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()


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
beta_x = rbc["beta_x"]
beta_y = rbc["beta_y"]

xMin = 0.0
xMax = 8.6
yMin = 0.0
yMax = 65.0

# Plotting
plt.figure(figsize=(10, 6))
plt.xscale("linear")
plt.yscale("linear")
plt.xlim([xMin, xMax])
# plt.ylim([yMin, yMax])
plt.xlabel("s (m)", fontsize=30)
plt.ylabel("CS Twiss beta (m)", fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.grid(True)

# Plot the data
plt.plot(s, beta_x, "b", label="Horizontal", linewidth=2, linestyle="solid")

plt.plot(s, beta_y, "r", label="Vertical", linewidth=2, linestyle="solid")

# Show plot
plt.legend(fontsize=20)

plt.tight_layout()
if args.save_png:
    plt.savefig("triplet.png")
else:
    plt.show()
