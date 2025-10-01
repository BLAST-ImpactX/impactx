#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import argparse

from plot_APL_ChrPlasmaLens import millimeter, plot_sigmas, plt, read_time_series

# options to run this script, this one is used by the CTest harness
parser = argparse.ArgumentParser(description="Plot the ChrPlasmaLens_zero benchmark.")
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

# import matplotlib.pyplot as plt

# read reduced diagnostics
rbc = read_time_series("diags/reduced_beam_characteristics.*")

# Plot beam transverse sizes
plot_sigmas(rbc)

# Start/end
plt.axhline(2.737665020201518e-05 * millimeter, ls="--", color="k")
# mid
plt.axhline(10e-6 * millimeter, ls="--", color="k")
plt.axvline(10e-3, ls="--", color="k")
if args.save_png:
    plt.savefig("APL_ChrPlasmaLens_zero-sigma.png")
else:
    plt.show()
