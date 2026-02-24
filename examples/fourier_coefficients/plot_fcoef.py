#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import argparse

import matplotlib.pyplot as plt
from fcoef import read_data
from matplotlib.ticker import MaxNLocator

# options to run this script
parser = argparse.ArgumentParser(
    description="Plot the reconstruction of on-axis data from Fourier coefficients."
)
parser.add_argument(
    "--save-png", action="store_true", help="non-interactive run: save to PNGs"
)
args = parser.parse_args()

data = read_data("onaxis_data.in")
z = data[:, 0]
g_onaxis = data[:, 1]
data = read_data("onaxis_data.out")
g_onaxis_reconstructed = data[:, 1]
gp_onaxis_reconstructed = data[:, 2]
gpp_onaxis_reconstructed = data[:, 3]
diff_onaxis = g_onaxis - g_onaxis_reconstructed

# print beam transverse size over steps
f = plt.figure(figsize=(9, 4.8))
ax1 = f.gca()
input_data = ax1.plot(z, g_onaxis, label=r"input data $g(z)$")
output_data = ax1.plot(z, g_onaxis_reconstructed, label=r"reconstructed data")
diff_data = ax1.plot(z, diff_onaxis, label=r"difference")

ax1.legend(
    handles=input_data + output_data + diff_data, loc="upper center", prop={"size": 11}
)
ax1.set_xlabel(r"$z$ [m]", fontsize=14)
ax1.set_ylabel(r"gradient", fontsize=14)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("fodo_sigma.png")
else:
    plt.show()

# print first derivative from reconstructed data
f = plt.figure(figsize=(9, 4.8))
ax1 = f.gca()
derivative1 = ax1.plot(z, gp_onaxis_reconstructed, label=r"$g'(z)$")

ax1.legend(handles=derivative1, loc="upper center", prop={"size": 11})
ax1.set_xlabel(r"$z$ [m]", fontsize=14)
ax1.set_ylabel(r"derivative", fontsize=14)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("fodo_polarization.png")
else:
    plt.show()

# print first and second derivatives from reconstructed data
f = plt.figure(figsize=(9, 4.8))
ax1 = f.gca()
derivative2 = ax1.plot(z, gpp_onaxis_reconstructed, label=r"$g''(z)$")

ax1.legend(handles=derivative2, loc="upper center", prop={"size": 11})
ax1.set_xlabel(r"$z$ [m]", fontsize=14)
ax1.set_ylabel(r"derivative", fontsize=14)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.tight_layout()
if args.save_png:
    plt.savefig("fodo_polarization.png")
else:
    plt.show()
