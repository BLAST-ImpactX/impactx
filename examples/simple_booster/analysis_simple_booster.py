#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
import math

import numpy as np
import openpmd_api as io
import scipy
from booster_impactx_lattice import get_lattice

from impactx import elements

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

# initial/final beam (TODO)
# series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
# last_step = list(series.iterations)[-1]
# initial = series.iterations[1].particles["beam"].to_df()
# final = series.iterations[last_step].particles["beam"].to_df()

# compare number of particles (TODO)
# num_particles = 10000
# assert num_particles == len(initial)
# assert num_particles == len(final)

# TODO for Eric: add more tests as needed, e.g., checking beam moments


series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

beta_ref = series.iterations[1].particles["beam"].get_attribute("beta_ref")

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

print("\t\tInitial\t\tFinal")
print("\t\t=======\t\t=====")
print(f"std x\t\t{sigma_x[0]:7g}\t{sigma_x[-1]:7g}")
print(f"std px\t\t{sigma_px[0]:7g}\t{sigma_px[-1]:7g}")
print(f"emittance x\t{emittance_x[0]:7g}\t{emittance_x[-1]:7g}")
print()
print(f"std y\t\t{sigma_y[0]:7g}\t{sigma_y[-1]:7g}")
print(f"std py\t\t{sigma_py[0]:7g}\t{sigma_py[-1]:7g}")
print(f"emittance y\t{emittance_y[0]:7g}\t{emittance_y[-1]:7g}")
print()
print(f"std t\t\t{sigma_t[0]:7g}\t\t{sigma_t[-1]:7g}")
print(f"std pt\t\t{sigma_pt[0]:7g}\t{sigma_pt[-1]:7g}")
print(f"emittance t\t{emittance_t[0]:7g}\t{emittance_t[-1]:7g}")

assert math.isclose(sigma_x[0], 0.0083783, rel_tol=1.0e-4)
assert math.isclose(sigma_x[-1], 0.0083783, rel_tol=1.0e-2), (
    "sigma_x change over one turn too large"
)

assert math.isclose(sigma_px[0], 0.000226298, rel_tol=1.0e-4)
assert math.isclose(sigma_px[-1], 0.000226298, rel_tol=1.0e-2), (
    "sigma_px change over one turn too large"
)

assert math.isclose(emittance_x[0], 1.89474e-06, rel_tol=1.0e-4)
assert math.isclose(emittance_x[1], 1.89474e-06, rel_tol=1.0e-2), (
    "emittance_x change over one turn too large"
)

assert math.isclose(sigma_y[0], 0.00299246, rel_tol=1.0e-4)
assert math.isclose(sigma_y[-1], 0.00299246, rel_tol=1.0e-2), (
    "sigma_y change over one turn too large"
)

assert math.isclose(sigma_py[0], 0.000575451, rel_tol=1.0e-4)
assert math.isclose(sigma_py[-1], 0.000575451, rel_tol=1.0e-2), (
    "sigma_py change over one turn too large"
)

assert math.isclose(emittance_y[0], 1.72198e-06, rel_tol=1.0e-4)
assert math.isclose(emittance_y[1], 1.72198e-06, rel_tol=1.0e-2), (
    "emittance_y change over one turn too large"
)

assert math.isclose(sigma_t[0], 1.16878, rel_tol=1.0e-4)
assert math.isclose(sigma_t[-1], 1.16878, rel_tol=1.0e-2), (
    "sigma_t change over one turn too large"
)

assert math.isclose(sigma_pt[0], 0.000914371, rel_tol=1.0e-4)
assert math.isclose(sigma_pt[-1], 0.000914371, rel_tol=1.0e-2), (
    "sigma_pt change over one turn too large"
)

assert math.isclose(emittance_t[0], 0.00106837, rel_tol=1.0e-4)
assert math.isclose(emittance_t[1], 0.00106837, rel_tol=1.0e-2), (
    "emittance_t change over one turn too large"
)
