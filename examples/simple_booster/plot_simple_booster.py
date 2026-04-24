#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import openpmd_api as io
import scipy

from booster_impactx_lattice import get_lattice

from impactx import ImpactX, elements

series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

beta_ref = series.iterations[1].particles["beam"].get_attribute("beta_ref")

# columns in rbc file
# step s mean_x min_x max_x mean_y min_y max_y mean_t min_t max_t sigma_x sigma_y sigma_t mean_px min_px max_px mean_py min_py max_py mean_pt min_pt max_pt sigma_px sigma_py sigma_pt emittance_x emittance_y emittance_t alpha_x alpha_y alpha_t beta_x beta_y beta_t dispersion_x dispersion_px dispersion_y dispersion_py emittance_xn emittance_yn emittance_tn charge_C

rbc = np.loadtxt("diags/reduced_beam_characteristics.0.0", skiprows=4)

step = rbc[:, 0]
s = rbc[:, 1]

mean_x = rbc[:, 2]
min_x = rbc[:, 3]
max_x = rbc[:, 4]

mean_y = rbc[:, 5]
min_y = rbc[:, 6]
max_y = rbc[:, 7]

mean_t = rbc[:, 8]
min_t = rbc[:, 9]
max_t = rbc[:, 10]

sigma_x = rbc[:, 11]
sigma_y = rbc[:, 12]
sigma_t = rbc[:, 13]

mean_px = rbc[:, 14]
min_px = rbc[:, 15]
max_px = rbc[:, 16]

mean_py = rbc[:, 17]
min_py = rbc[:, 18]
max_py = rbc[:, 19]

mean_pt = rbc[:, 20]
min_pt = rbc[:, 21]
max_pt = rbc[:, 22]

sigma_px = rbc[:, 23]
sigma_py = rbc[:, 24]
sigma_pt = rbc[:, 25]

emittance_x = rbc[:, 26]
emittance_y = rbc[:, 27]
emittance_t = rbc[:, 28]

alpha_x = rbc[:, 29]
alpha_y = rbc[:, 30]
alpha_t = rbc[:, 31]

beta_x = rbc[:, 32]
beta_y = rbc[:, 33]
beta_t = rbc[:, 34]

dispersion_x = rbc[:, 35]
dispersion_px = rbc[:, 36]
dispersion_y = rbc[:, 37]
dispersion_py = rbc[:, 38]

emittance_xn = rbc[:, 39]
emittance_yn = rbc[:, 40]
emittance_tn = rbc[:, 41]

charge = rbc[:, 42] / scipy.constants.eV

f, ax = plt.subplots(2, 2)
plt.suptitle("sigmas vs. s")
ax[0, 0].plot(s, sigma_x*1000, label="sigma_x [mm]")
ax[0, 0].legend(loc="best")
ax[0, 1].plot(s, sigma_y*1000, label="sigma_y [mm]")
ax[0, 1].legend(loc="best")
ax[1, 0].plot(s, sigma_t, label="sigma_t [m]")
ax[1, 0].legend(loc="best")
ax[1, 1].plot(s, sigma_pt*1000, label="sigma_pt [x1000]")
ax[1, 1].legend(loc='best')

plt.figure()
plt.title("charge")
plt.plot(s, charge, label="charge [C]")
plt.legend(loc="best")

f, ax = plt.subplots(3, 2)
ax[0, 0].plot(s, min_x*1000, label="min x [mm]")
ax[0, 0].legend(loc="best")

ax[0, 1].plot(s, max_x*1000, label="max x [mm]")
ax[0, 1].legend(loc="best")

ax[1, 0].plot(s, min_y*1000, label="min y [mm]")
ax[1, 0].legend(loc="best")

ax[1, 1].plot(s, max_y*1000, label="max y [mm]")
ax[1, 1].legend(loc="best")

ax[2, 0].plot(s, min_t*beta_ref, label="min z [m]")
ax[2, 0].legend(loc="best")

ax[2, 1].plot(s, max_t*beta_ref, label="max z [m]")
ax[2, 1].legend(loc="best")


f, ax = plt.subplots(3, 2)

ax[0, 0].plot(initial["position_x"]*1000, initial["momentum_x"]*1000, '.', label="initial x [mm] vs. px [mr]")
ax[0, 0].legend(loc="best")

ax[0, 1].plot(final["position_x"]*1000, final["momentum_x"]*1000, '.', label="final x [mm] vs. px [mr]")
ax[0, 1].legend(loc="best")

ax[1, 0].plot(initial["position_y"]*1000, initial["momentum_y"]*1000, '.', label="initial y [mm] vs. py [mr]")
ax[1, 0].legend(loc="best")

ax[1, 1].plot(final["position_y"]*1000, final["momentum_y"]*1000, '.', label="final y [mm] vs. py [mr]")
ax[1, 1].legend(loc="best")

ax[2, 0].plot(initial["position_t"], initial["momentum_t"]*1000, '.', label="initial t [m] vs. pt [x1000]")
ax[2, 0].legend(loc="best")

ax[2, 1].plot(final["position_t"], final["momentum_t"]*1000, '.', label="final t [m] vs. pt [x1000]")
ax[2, 1].legend(loc="best")

plt.show()
