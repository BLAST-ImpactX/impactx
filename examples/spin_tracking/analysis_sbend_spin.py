#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io


def get_polarization(beam):
    """Calculate polarization vector, given by the mean values of spin components.

    Returns
    -------
    polarization_x, polarization_y, polarization_z
    """
    polarization_x = np.mean(beam["spin_x"])
    polarization_y = np.mean(beam["spin_y"])
    polarization_z = np.mean(beam["spin_z"])

    return (polarization_x, polarization_y, polarization_z)


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

# compare number of particles
num_particles = 100000
assert num_particles == len(initial)
assert num_particles == len(final)

# numerical parameters based on input file

gyro_anomaly = 0.001159652181644  # for electrons
quad_gradient = 100  # value in 1/m^2 from input
sigma_x = 0.003  # value in m = lambdaX from input
sigma_px = 0.2  # value in rad = lambdaPx from input
sigma_pt = 0.2  # value = lambdaPt from input
Pxi = 0.4  # polarization_x from input
Pyi = 0.9  # polarization_y from input
Pzi = 0.1  # polarization_z from input

rel_gamma = 3914.902367118367246  # relativistic gamma
rel_beta = 0.999999967376681  # relativistic beta
h = 1.0 / 10.0  # 1/radius of bend curvature from input

print("Initial Beam:")
polarization_x, polarization_y, polarization_z = get_polarization(initial)
print(
    f"  polarization_x={polarization_x:e} polarization_y={polarization_y:e} polarization_z={polarization_z:e}"
)

atol = 1.3 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  atol={atol}")

assert np.allclose(
    [polarization_x, polarization_y, polarization_z],
    [
        Pxi,
        Pyi,
        Pzi,
    ],
    atol=atol,
)

# predicted final polarization
gyro_const = 1.0 + gyro_anomaly * rel_gamma
sin1 = np.sin(np.pi / (gyro_anomaly * rel_gamma))
sin2 = np.sin(2 * np.pi / (gyro_anomaly * rel_gamma))
damping_eigenvalue2 = (
    4 * sigma_px**2 * gyro_const**2 * sin1**4
    + h**2 * gyro_const**2 * sigma_x**2 * sin2**2
    + sigma_pt**2 * (-2.0 * np.pi * rel_beta**2 + gyro_const * sin2) ** 2 / rel_beta**2
)
damping_eigenvalue = np.sqrt(damping_eigenvalue2)
damping_factor = np.exp(-(damping_eigenvalue**2) / 2.0)

Pxf = damping_factor * Pxi
Pyf = Pyi
Pzf = damping_factor * Pzi

print("")
print("Predicted Final Polarization:")
print(f"  polarization_x={Pxf:e} polarization_y={Pyf:e} polarization_z={Pzf:e}")

print("")
print("Final Beam:")
polarization_x, polarization_y, polarization_z = get_polarization(final)
print(
    f"  polarization_x={polarization_x:e} polarization_y={polarization_y:e} polarization_z={polarization_z:e}"
)

atol = 1.3 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  atol={atol}")

assert np.allclose(
    [polarization_x, polarization_y, polarization_z],
    [
        Pxf,
        Pyf,
        Pzf,
    ],
    atol=atol,
)
