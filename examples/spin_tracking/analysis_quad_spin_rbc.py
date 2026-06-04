#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd


# Load data from envelope simulation
def read_time_series(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.

    Returns
    -------
    pandas.DataFrame
    """

    def read_file(file_pattern):
        for filename in glob.glob(file_pattern):
            df = pd.read_csv(filename, delimiter=r"\s+")
            if "step" not in df.columns:
                step = int(re.findall(r"[0-9]+", filename)[0])
                df["step"] = step
            yield df

    return pd.concat(
        read_file(file_pattern),
        axis=0,
        ignore_index=True,
    )  # .set_index('id')


def data_is_double(file_pattern):
    """Detect float precision (single vs double) from a text diagnostic's digits."""
    text = "".join(Path(f).read_text() for f in glob.glob(file_pattern))
    return any(
        len(m.replace(".", "") if e else m.replace(".", "").strip("0")) >= 12
        for m, e in re.findall(r"(\d*\.\d+)([eE][+-]?\d+)?", text)
    )


rbc = read_time_series("diags/reduced_beam_characteristics.*")
is_double = data_is_double("diags/reduced_beam_characteristics.*")

# numerical parameters based on input file
gryo_anomaly = 0.001159652181644  # for electrons
rel_gamma = 196.69511809100055  # for 100 MeV
quad_gradient = 100  # value in 1/m^2 from input
sigma_y = 0.0003  # value in m = lambdaY from input
sigma_py = 0.0002  # value in rad = lambdaPy from input
Pxi = 0.4  # polarization_x from input
Pyi = 0.9  # polarization_y from input
Pzi = 0.1  # polarization_z from input

print("Initial Beam:")
polarization_x = rbc["mean_sx"].iloc[0]
polarization_y = rbc["mean_sy"].iloc[0]
polarization_z = rbc["mean_sz"].iloc[0]
print(
    f"  polarization_x={polarization_x:e} polarization_y={polarization_y:e} polarization_z={polarization_z:e}"
)

# Number of particles (from input file)
num_particles = 100000

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
polarization_x = rbc["mean_sx"].iloc[-1]
polarization_y = rbc["mean_sy"].iloc[-1]
polarization_z = rbc["mean_sz"].iloc[-1]
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

# numerical tests of spin moment conditions
sigma_sx = rbc["sigma_sx"].iloc[-1]
sigma_sy = rbc["sigma_sy"].iloc[-1]
sigma_sz = rbc["sigma_sz"].iloc[-1]
polarization = np.sqrt(polarization_x**2 + polarization_y**2 + polarization_z**2)
condition = sigma_sx**2 + sigma_sy**2 + sigma_sz**2 + polarization**2

print("")
print(f"Spin moment consistency condition = {condition:e}")

atol = 1.0e-12 if is_double else 1.0e-3
print(f"  atol={atol}")

assert np.allclose(
    [condition],
    [
        1.0,
    ],
    atol=atol,
)
