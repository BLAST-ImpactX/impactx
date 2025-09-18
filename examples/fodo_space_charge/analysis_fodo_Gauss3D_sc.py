#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Ji Qiang
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

print("Initial Beam:")
sig_xi, sig_yi, sig_ti, emittance_xi, emittance_yi, emittance_ti = get_moments(initial)
print(f"  sigx={sig_xi:e} sigy={sig_yi:e} sigt={sig_ti:e}")
print(
    f"  emittance_x={emittance_xi:e} emittance_y={emittance_yi:e} emittance_t={emittance_ti:e}"
)

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sig_xi, sig_yi, sig_ti, emittance_xi, emittance_yi, emittance_ti],
    [7.515765e-05, 7.511883e-05, 9.997395e-4, 2.001510e-09, 1.999755e-09, 1.999289e-06],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam:")
sig_xf, sig_yf, sig_tf, emittance_xf, emittance_yf, emittance_tf = get_moments(initial)
print(f"  sigx={sig_xf:e} sigy={sig_yf:e} sigt={sig_tf:e}")
print(
    f"  emittance_x={emittance_xf:e} emittance_y={emittance_yf:e} emittance_t={emittance_tf:e}"
)

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sig_xf, sig_yf, sig_tf, emittance_xf, emittance_yf, emittance_tf],
    [
        7.51576586332169e-05,
        7.511883208451813e-05,
        0.0009997395499750136,
        2.0015106608723994e-09,
        1.999755254276969e-09,
        1.9992898444562777e-06,
    ],
    rtol=rtol,
    atol=atol,
)
