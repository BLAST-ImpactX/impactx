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
import openpmd_api as io
import pandas as pd
from matplotlib.ticker import MaxNLocator


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


# options to run this script
parser = argparse.ArgumentParser(description="Plot action of the polygon aperture.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()


f, axs = plt.subplots(1, 2)

axs[0].scatter(initial['position_x']*1.0e3,
               initial['position_y']*1.0e3)
axs[0].set_title("initial")
axs[0].set_xlabel(r"$x$ [mm]")
axs[0].set_ylabel(r"$y$ [mm]")
axs[0].set_xlim([-5.5, 5.5])
axs[0].set_ylim([-5.5, 5.5])

axs[1].scatter(final['position_x']*1.0e3,
                final['position_y']*1.0e3)
axs[1].set_title("final")
axs[1].set_xlabel(r"$x$ [mm]")
axs[1].set_ylabel(r"$y$ [mm]")
axs[1].set_xlim([-5.5, 5.5])
axs[1].set_ylim([-5.3, 5.3])
                

plt.tight_layout()
if args.save_png:
    plt.savefig("polygon_aperture.png")
else:
    plt.show()
