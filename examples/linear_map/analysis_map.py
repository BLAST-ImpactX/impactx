#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
from scipy.stats import moment


def get_moments(beam):
    """Calculate standard deviations of beam position & momenta
    and emittance values

    Returns
    -------
    sigx, sigy, sigt, emittance_x, emittance_y, emittance_t
    """
    sigx = moment(beam["position_x"], moment=2) ** 0.5  # variance -> std dev.
    sigpx = moment(beam["momentum_x"], moment=2) ** 0.5
    sigy = moment(beam["position_y"], moment=2) ** 0.5
    sigpy = moment(beam["momentum_y"], moment=2) ** 0.5
    sigt = moment(beam["position_t"], moment=2) ** 0.5
    sigpt = moment(beam["momentum_t"], moment=2) ** 0.5

    epstrms = beam.cov(ddof=0)
    emittance_x = (sigx**2 * sigpx**2 - epstrms["position_x"]["momentum_x"] ** 2) ** 0.5
    emittance_y = (sigy**2 * sigpy**2 - epstrms["position_y"]["momentum_y"] ** 2) ** 0.5
    emittance_t = (sigt**2 * sigpt**2 - epstrms["position_t"]["momentum_t"] ** 2) ** 0.5

    return (sigx, sigy, sigt, emittance_x, emittance_y, emittance_t)


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
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(initial)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        6.363961030678928e-6,
        28.284271247461902e-9,
        0.0035,
        0.27e-9,
        1.0e-12,
        1.33e-6,
    ],
    rtol=rtol,
    atol=atol,
)

print("")
print("Final Beam:")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(final)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        6.363961030678928e-6,
        28.284271247461902e-9,
        0.0035,
        0.27e-9,
        1.0e-12,
        1.33e-6,
    ],
    rtol=rtol,
    atol=atol,
)

# Specify time series for particle j
j = 5
print(f"output for particle index = {j}")

# Create array of TBT data values
x = []
px = []
y = []
py = []
t = []
pt = []
n = 0
for k_i, i in series.iterations.items():
    beam = i.particles["beam"]
    turn = beam.to_df()
    x.append(turn["position_x"][j])
    px.append(turn["momentum_x"][j])
    y.append(turn["position_y"][j])
    py.append(turn["momentum_y"][j])
    t.append(turn["position_t"][j])
    pt.append(turn["momentum_t"][j])
    n = n + 1

# Output number of periods in data series
nturns = len(x)
print(f"number of periods = {nturns}")
print()

# Approximate the tune and closed orbit using the 4-turn formula:

# from x data only
argument = (x[0] - x[1] + x[2] - x[3]) / (2.0 * (x[1] - x[2]))
tunex = np.arccos(argument) / (2.0 * np.pi)
print(f"tune output from 4-turn formula, using x data = {tunex}")

# from y data only
argument = (y[0] - y[1] + y[2] - y[3]) / (2.0 * (y[1] - y[2]))
tuney = np.arccos(argument) / (2.0 * np.pi)
print(f"tune output from 4-turn formula, using y data = {tuney}")

# from t data only
argument = (t[0] - t[1] + t[2] - t[3]) / (2.0 * (t[1] - t[2]))
tunet = np.arccos(argument) / (2.0 * np.pi)
print(f"tune output from 4-turn formula, using t data = {tunet}")

rtol = 1.0e-3
print(f"  rtol={rtol}")

assert np.allclose(
    [tunex, tuney, tunet],
    [
        0.139,
        0.219,
        0.0250,
    ],
    rtol=rtol,
)
