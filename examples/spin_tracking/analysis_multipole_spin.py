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
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

print("Compare initial and final beam moments:")

sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(initial)
sigxf, sigyf, sigtf, emittance_xf, emittance_yf, emittance_tf = get_moments(final)

atol = 0.0  # ignored
rtol = 3.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [sigxf, sigyf, sigtf, emittance_xf, emittance_yf, emittance_tf],
    rtol=rtol,
    atol=atol,
)

print("Compare initial and final spins:")

sxi = initial["spin_x"]
syi = initial["spin_y"]
szi = initial["spin_z"]

sxf = final["spin_x"]
syf = final["spin_y"]
szf = final["spin_z"]

dspin2 = (sxf - sxi) ** 2 + (syf - syi) ** 2 + (szf - szi) ** 2
dspin = np.sqrt(dspin2)
dspinmax = dspin.max()

print("Change in the spin:")
print("||delta s||_max", dspinmax)

atol = 5.0e-8
print(f"  atol={atol}")

assert np.allclose(
    [dspinmax],
    [
        0.0,
    ],
    atol=atol,
)
