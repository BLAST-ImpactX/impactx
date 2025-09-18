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
series1 = io.Series("diags/openPMD/TCPhosphor.h5", io.Access.read_only)
last_step = list(series1.iterations)[-1]
beam_initial = series1.iterations[last_step].particles["beam"]
initial = beam_initial.to_df()

series2 = io.Series("diags/openPMD/UC_VisaEBeam8.h5", io.Access.read_only)
last_step = list(series2.iterations)[-1]
beam_final = series2.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

gamma_ref = beam_initial.get_attribute("gamma_ref")
bg = np.sqrt(gamma_ref**2 - 1.0)
print("Initial Beam (at TCPhosphor screen):")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(initial)

emittance_xn = emittance_x * bg
emittance_yn = emittance_y * bg
emittance_tn = emittance_t * bg

print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_xn={emittance_xn:e} emittance_yn={emittance_yn:e} emittance_tn={emittance_tn:e}"
)

atol = 0.0  # ignored
rtol = 20.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_xn, emittance_yn, emittance_tn],
    [
        4.786098e-04,
        2.983777e-04,
        1.147205e-06,
        3.791183e-06,
        2.164368e-06,
        5.589590e-06,
    ],
    rtol=rtol,
    atol=atol,
)

gamma_ref = beam_final.get_attribute("gamma_ref")
bg = np.sqrt(gamma_ref**2 - 1.0)

print("")
print("Final Beam (at VisaEbeam8 screen):")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(final)

emittance_xn = emittance_x * bg
emittance_yn = emittance_y * bg
emittance_tn = emittance_t * bg

print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_xn={emittance_xn:e} emittance_yn={emittance_yn:e} emittance_tn={emittance_tn:e}\n"
)

atol = 0.0  # ignored
rtol = 20.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_xn, emittance_yn, emittance_tn],
    [
        7.554152e-03,
        2.151019e-03,
        9.054048e-04,
        6.090367e-03,
        6.238249e-04,
        4.425367e-03,
    ],
    rtol=rtol,
    atol=atol,
)
