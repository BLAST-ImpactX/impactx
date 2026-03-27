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
gryo_anomaly = 0.001159652181644  # for electrons
rel_gamma = beam_initial.get_attribute("gamma_ref")  # for 100 MeV
quad_gradient = 100  # value in 1/m^2 from input
sigma_y = 0.0003  # value in m = lambdaY from input
sigma_py = 0.0002  # value in rad = lambdaPy from input
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
damping_eigenvalue = (1 + gryo_anomaly * rel_gamma) * np.sqrt(
    sigma_py**2 * (np.cosh(2 * np.pi) - 1) ** 2
    + sigma_y**2 * quad_gradient * np.sinh(2 * np.pi) ** 2
)
damping_factor = np.exp(-(damping_eigenvalue**2) / 2.0)
Pxf = Pxi
Pyf = damping_factor * Pyi
Pzf = damping_factor * Pzi

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
