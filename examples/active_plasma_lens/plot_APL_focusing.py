#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import argparse

from analysis_APL import read_time_series
from plot_APL import millimeter, plot_sigmas, plt
from run_APL import analytic_sigma_function

# options to run this script, this one is used by the CTest harness
parser = argparse.ArgumentParser(
    description="Plot the ChrPlasmaLens_focusing and ConstF_tracking_focusing benchmarks."
)
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

# import matplotlib.pyplot as plt

# read reduced diagnostics
rbc = read_time_series("diags/reduced_beam_characteristics.*")

# Plot beam transverse sizes
plot_sigmas(rbc)

# Analytical estimates
# Start/end
plt.axhline(7.161196476484095e-05 * millimeter, ls="--", color="k")
plt.axhline(100e-6 * millimeter, ls="--", color="k")
# plt.axvline(10e-3, ls="--", color="k")
# As function of s
(s, sigma) = analytic_sigma_function(-1000, 100e-6)
plt.plot(s, sigma * millimeter, ls="--", color="green", label="Analytical")

plt.legend(loc="center left")
plt.title(r"Focusing e$^-$, 200 MeV, $g$ = -1000 [T/m]")
plt.tight_layout()

if args.save_png:
    plt.savefig("APL_focusing-sigma.png")
else:
    plt.show()
