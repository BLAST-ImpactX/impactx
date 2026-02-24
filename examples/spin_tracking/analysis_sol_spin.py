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
beam_initial = series.iterations[1].particles["beam"]
gyro_anomaly = 1.7928473446  # for protons
rel_beta = beam_initial.get_attribute("beta_ref")
sol_strength_ks = 1  # value in 1/m^2 from input
sol_length_ds = 3.504584663108292
sigma_pt = 0.2  # value in rad = lambdaPt from input
Pxi = 0.4  # polarization_x from input
Pyi = 0.9  # polarization_y from input
Pzi = 0.1  # polarization_z from input

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
damping_eigenvalue = (
    (1 + gyro_anomaly) * sol_strength_ks * sol_length_ds * sigma_pt / rel_beta
)
damping_factor = np.exp(-(damping_eigenvalue**2) / 2.0)
cosG = np.cos(2 * (1 + gyro_anomaly) * np.pi / gyro_anomaly)
sinG = np.sin(2 * (1 + gyro_anomaly) * np.pi / gyro_anomaly)
Pxf = damping_factor * (Pxi * cosG + Pyi * sinG)
Pyf = damping_factor * (-Pxi * sinG + Pyi * cosG)
Pzf = Pzi

print("")
print("Predicted Final Polarization:")
print(f"  polarization_x={Pxf:e} polarization_y={Pyf:e} polarization_z={Pzf:e}")

print("")
print("Final Beam:")
polarization_x, polarization_y, polarization_z = get_polarization(final)
print(
    f"  polarization_x={polarization_x:e} polarization_y={polarization_y:e} polarization_z={polarization_z:e}"
)

atol = 2.0 * num_particles**-0.5  # from random sampling of a smooth distribution
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
