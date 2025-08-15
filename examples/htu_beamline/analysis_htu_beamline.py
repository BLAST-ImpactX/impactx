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
initial = series1.iterations[last_step].particles["beam"].to_df()

series2 = io.Series("diags/openPMD/UC_VisaEBeam8.h5", io.Access.read_only)
last_step = list(series2.iterations)[-1]
beam_final = series2.iterations[last_step].particles["beam"]
final = beam_final.to_df()

# compare number of particles
num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

print("Initial Beam (at TCPhosphor screen):")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(initial)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 10.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        3.536004e-04,
        2.152397e-04,
        1.530077e-06,
        3.734497e-08,
        1.506135e-08,
        3.729656e-08,
    ],
    rtol=rtol,
    atol=atol,
)


print("")
print("Final Beam (at VisaEbeam8 screen):")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(final)
s_ref = beam_final.get_attribute("s_ref")
gamma_ref = beam_final.get_attribute("gamma_ref")
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}\n"
    f"  s_ref={s_ref:e} gamma_ref={gamma_ref:e}"
)

atol = 0.0  # ignored
rtol = 10.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, s_ref, gamma_ref],
    [
        1.759036e-03,
        2.337416e-04,
        1.175086e-04,
        2.813148e-06,
        6.167234e-07,
        2.858550e-06,
        9.459980,
        1.956951e02,
    ],
    rtol=rtol,
    atol=atol,
)
