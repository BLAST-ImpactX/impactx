#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import openpmd_api as io
from scipy.stats import moment


def get_moments(beam):
    """Calculate standard deviations of beam position & momenta
    and emittance values

    Returns
    -------
    sigx, sigy, sigt, sigpx, sigpy, sigpt, emittance_x, emittance_y, emittance_t
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

    return (
        sigx,
        sigy,
        sigt,
        sigpx,
        sigpy,
        sigpt,
        emittance_x,
        emittance_y,
        emittance_t,
    )


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial_beam = series.iterations[1].particles["beam"]
final_beam = series.iterations[last_step].particles["beam"]
initial = initial_beam.to_df()
final = final_beam.to_df()

# compare number of particles
num_particles = 100000
assert num_particles == len(initial)

# problem parameters
beam_radius = 2.0e-3
aperture_radius = 3.5e-3
correlation_k = 0.5
drift_distance = 6.0

print("Initial Beam:")
sigx, sigy, sigt, sigpx, sigpy, sigpt, emittance_x, emittance_y, emittance_t = (
    get_moments(initial)
)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(f"  sigpx={sigpx:e} sigpy={sigpy:e} sigpt={sigpt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 2.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt],
    [
        beam_radius / 2.0,
        beam_radius / 2.0,
        beam_radius / 2.0,
    ],
    rtol=rtol,
    atol=atol,
)

atol = 1.0e-11  # ignored
print(f"  atol={atol}")

assert np.allclose(
    [sigpx, sigpy, sigpt, emittance_x, emittance_y, emittance_t],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    atol=atol,
)


# calculation of predicted final beam parameters
beam_radius_no_aperture = beam_radius * (1.0 + correlation_k * drift_distance)
beam_radius_with_aperture = min(beam_radius_no_aperture, aperture_radius)

fractional_loss = 1.0 - min(1.0, (aperture_radius / beam_radius_no_aperture) ** 2)
sigma_x_final = beam_radius_with_aperture / 2.0
sigma_px_final = correlation_k / (1.0 + correlation_k * drift_distance) * sigma_x_final

print("")
print("Predicted Final Beam:")
print(f"  sigx={sigma_x_final:e} sigy={sigma_x_final:e} sigt={beam_radius / 2.0:e}")
print(f"  sigpx={sigma_px_final:e} sigpy={sigma_px_final:e} sigpt=0.0")
print(f"  fractional_loss={fractional_loss:e}")


print("")
print("Final Beam:")
sigx, sigy, sigt, sigpx, sigpy, sigpt, emittance_x, emittance_y, emittance_t = (
    get_moments(final)
)
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(f"  sigpx={sigpx:e} sigpy={sigpy:e} sigpt={sigpt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

atol = 0.0  # ignored
rtol = 2.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sigx, sigy, sigt, sigpx, sigpy],
    [
        sigma_x_final,
        sigma_x_final,
        beam_radius / 2.0,
        sigma_px_final,
        sigma_px_final,
    ],
    rtol=rtol,
    atol=atol,
)

charge_i = initial_beam.get_attribute("charge_C")
charge_f = final_beam.get_attribute("charge_C")

loss_pct = 100.0 * (charge_i - charge_f) / charge_i

print(f" fractional loss (%) = {loss_pct}")

atol = 0.2  # tolerance 0.2%
print(f"  atol={atol}")
assert np.allclose(
    [loss_pct],
    [100 * fractional_loss],
    atol=atol,
)
