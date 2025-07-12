#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import numpy as np
import openpmd_api as io
from scipy.constants import c, e, m_e
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

initial_step = list(series.iterations)[0]
last_step = list(series.iterations)[-1]
initial = series.iterations[initial_step].particles["beam"].to_df()
beam_final = series.iterations[last_step].particles["beam"]
final = beam_final.to_df()

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
rtol = 3.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        1.768569e-04,
        1.195967e-04,
        1.034426e-06,
        1.356459e-08,
        9.374827e-09,
        2.581087e-08,
    ],
    rtol=rtol,
    atol=atol,
)

print("")
print("Final Beam:")
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(final)
s_ref = beam_final.get_attribute("s_ref")
gamma_ref = beam_final.get_attribute("gamma_ref")
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}\n"
    f"  s_ref={s_ref:e} gamma_ref={gamma_ref:e}"
)

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        2.458574e-04,
        1.513417e-04,
        1.068497e-06,
        1.797301e-08,
        1.009999e-08,
        2.661524e-08,
    ],
    rtol=rtol,
    atol=atol,
)


# join tables on particle ID, so we can compare the same particle initial->final
beam_joined = final.join(initial, lsuffix="_final", rsuffix="_initial")
xi = beam_joined["position_x_initial"]
yi = beam_joined["position_y_initial"]
pxi = beam_joined["momentum_x_initial"]
pyi = beam_joined["momentum_y_initial"]
pti = beam_joined["momentum_t_initial"]
xf = beam_joined["position_x_final"]
yf = beam_joined["position_y_final"]
pxf = beam_joined["momentum_x_final"]
pyf = beam_joined["momentum_y_final"]
ptf = beam_joined["momentum_t_final"]

# Parameters appearing in the Hamiltonian
beta_ref = beam_final.get_attribute("beta_ref")
bg_ref = beam_final.get_attribute("beta_gamma_ref")
rigidity = m_e * bg_ref * c / e
B_gradient = 207.0
k = B_gradient / rigidity

# Evaluate the change in Hamiltonian value for each particle
Hi = (
    -np.sqrt(1.0 - 2.0 / beta_ref * pti + pti**2 - pxi**2 - pyi**2)
    + (k / 2.0) * (xi**2 - yi**2)
    - pti / beta_ref
    + 1.0
)
Hf = (
    -np.sqrt(1.0 - 2.0 / beta_ref * ptf + ptf**2 - pxf**2 - pyf**2)
    + (k / 2.0) * (xf**2 - yf**2)
    - ptf / beta_ref
    + 1.0
)

H_sigma = (np.std(Hi) + np.std(Hf)) / 2.0
dH = (Hf - Hi).abs()
dH_max_relative = dH.max() / H_sigma

# particle-wise comparison of H & I initial to final
atol = 1.0e-3
rtol = 0.0  # large number
print()
print(f"  atol={atol} (ignored: rtol~={rtol})")

print(f"  Standard deviation in H: H_sigma = {H_sigma}")
print(f"  dH_max relative to H_sigma = {dH_max_relative}")
assert np.allclose(dH_max_relative, 0.0, rtol=rtol, atol=atol)
