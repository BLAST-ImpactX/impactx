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
    description="Plot the ChrPlasmaLens_zero and ConstF_tracking_zero benchmarks."
)
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

# read reduced diagnostics
rbc = read_time_series("diags/reduced_beam_characteristics.*")

# Plot beam transverse sizes
plot_sigmas(rbc)

# Analytical estimates
# Start/end
plt.axhline(2.737665020201518e-05 * millimeter, ls="--", color="k")
# mid
plt.axhline(10e-6 * millimeter, ls="--", color="k")
plt.axvline(10e-3, ls="--", color="k")
# As function of s
(s, sigma) = analytic_sigma_function(0.0, 10e-6)
plt.plot(s, sigma * millimeter, ls="--", color="green", label="Analytical")

plt.legend(loc="center")
plt.title(r"No-field e$^-$, 200 MeV, $g$ = 0 [T/m]")
plt.tight_layout()

if args.save_png:
    plt.savefig("APL_zero-sigma.png")
else:
    plt.show()
