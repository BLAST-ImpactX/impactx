#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import glob
import re

import matplotlib.pyplot as plt
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


# scaling to units
millimeter = 1.0e3  # m->mm
mrad = 1.0e3  # ImpactX uses "static units": momenta are normalized by the magnitude of the momentum of the reference particle p0: px/p0 (rad)
# mm_mrad = 1.e6
nm_rad = 1.0e9


def plot_sigmas(rbc):
    s = rbc["s"]
    sigma_x = rbc["sigma_x"] * millimeter
    sigma_y = rbc["sigma_y"] * millimeter

    # emittance_x = rbc["emittance_x"] * nm_rad
    # emittance_y = rbc["emittance_y"] * nm_rad

    # print beam transverse size over steps
    f = plt.figure()
    ax1 = f.gca()
    im_sigx = ax1.plot(s, sigma_x, label=r"$\sigma_x$")
    im_sigy = ax1.plot(s, sigma_y, label=r"$\sigma_y$")
    # ax2 = ax1.twinx()
    # ax2.set_prop_cycle(None)  # reset color cycle
    # im_emittance_x = ax2.plot(s, emittance_x, ":", label=r"$\epsilon_x$")
    # im_emittance_y = ax2.plot(s, emittance_y, ":", label=r"$\epsilon_y$")

    # ax1.legend(
    #    handles=im_sigx + im_sigy + im_emittance_x + im_emittance_y, loc="lower center"
    # )
    ax1.legend(handles=im_sigx + im_sigy, loc="best")
    ax1.set_xlabel(r"$z$ [m]")
    ax1.set_ylabel(r"$\sigma_{x,y}$ [mm]")
    # ax2.set_ylabel(r"$\epsilon_{x,y}$ [nm]")
    # ax2.set_ylim([1.5, 2.5])
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.tight_layout()

    # return (ax1,ax2)
