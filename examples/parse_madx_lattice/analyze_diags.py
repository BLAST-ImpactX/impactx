#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import openpmd_api as io
import scipy

# columns in rbc file
# step s mean_x min_x max_x mean_y min_y max_y mean_t min_t max_t sigma_x sigma_y sigma_t mean_px min_px max_px mean_py min_py max_py mean_pt min_pt max_pt sigma_px sigma_py sigma_pt emittance_x emittance_y emittance_t alpha_x alpha_y alpha_t beta_x beta_y beta_t dispersion_x dispersion_px dispersion_y dispersion_py emittance_xn emittance_yn emittance_tn charge_C


series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
iterations = list(series.iterations)
iter0 = iterations[0]
iter1 = iterations[-1]

s_by_iterations = np.array(
    [
        series.iterations[it].particles["beam"].get_attribute("s_ref")
        for it in iterations
    ]
)
npart_by_iterations = np.array(
    [
        series.iterations[it].particles["beam"]["position"]["x"].shape[0]
        for it in iterations
    ]
)

plt.figure()
plt.title("macroparticles")
plt.plot(s_by_iterations, npart_by_iterations)
plt.xlabel("s [m]")
plt.ylabel("macroparticles")


rbc_0 = np.loadtxt("diags/reduced_beam_characteristics.0.0", skiprows=1, max_rows=1)

rbc_1 = np.loadtxt("diags/reduced_beam_characteristics.0.0", skiprows=4)

rbc = np.vstack((rbc_0, rbc_1))

step_s = rbc[:, 0]
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


print("initial rbc stdx: ", sigma_x[0])
print(
    "initial monitor stdx: ",
    series.iterations[iter0].particles["beam"].get_attribute("sigma_x"),
)
print()
print("final rbc std: ", sigma_x[-1])
print(
    "final monitor stdx: ",
    series.iterations[iter1].particles["beam"].get_attribute("sigma_x"),
)

print("sigma_x")
print(sigma_x[:30])
print()
print("sigma_y")
print(sigma_y[:30])
print()

plt.figure()
plt.subplot(221)
plt.suptitle("sigmas vs. slice")
plt.plot(sigma_x, label="sigma_x")
plt.legend(loc="best")
plt.subplot(222)
plt.plot(sigma_y, label="sigma_y")
plt.legend(loc="best")
plt.subplot(223)
plt.plot(sigma_t, label="sigma_t")
plt.legend(loc="best")
plt.subplot(224)
plt.plot(charge, label="charge")
plt.legend(loc="best")

plt.figure()
plt.suptitle("sigmas vs. s")
plt.subplot(221)
plt.plot(s, sigma_x, label="sigma_x")
plt.legend(loc="best")
plt.subplot(222)
plt.plot(s, sigma_y, label="sigma_y")
plt.legend(loc="best")
plt.subplot(223)
plt.plot(s, sigma_t, label="sigma_t")
plt.legend(loc="best")
plt.subplot(224)
plt.plot(s, charge, label="charge")
plt.legend(loc="best")

plt.show()
