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
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

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
        0.7,
        0.0,
        0.0,
    ],
    atol=atol,
)

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
        0.7,
        0.0,
        0.0,
    ],
    atol=atol,
)
