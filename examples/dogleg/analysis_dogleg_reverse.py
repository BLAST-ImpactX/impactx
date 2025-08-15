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


def get_twiss(openpmd_beam):
    """Return Twiss functions from an openPMD particle species

    Returns
    -------
    alpha_x, beta_x, alpha_y, beta_y, dispersion_x, dispersion_px
    """
    alpha_x = openpmd_beam.get_attribute("alpha_x")
    beta_x = openpmd_beam.get_attribute("beta_x")
    d_x = openpmd_beam.get_attribute("dispersion_x")
    d_px = openpmd_beam.get_attribute("dispersion_px")
    alpha_y = openpmd_beam.get_attribute("alpha_y")
    beta_y = openpmd_beam.get_attribute("beta_y")
    # d_y = openpmd_beam.get_attribute("dispersion_y")
    # d_py = openpmd_beam.get_attribute(["dispersion_py")

    return (alpha_x, beta_x, alpha_y, beta_y, d_x, d_px)


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial_beam = series.iterations[1].particles["beam"]
initial = initial_beam.to_df()
final_beam = series.iterations[last_step].particles["beam"]
final = final_beam.to_df()

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
        1.924711e-03,
        2.165646e-05,
        1.102534e-04,
        7.668809e-09,
        1.018986e-10,
        8.588054e-09,
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
        2.057123e-05,
        6.911405e-05,
        2.012573e-05,
        8.174618e-11,
        1.018986e-10,
        1.151058e-08,
    ],
    rtol=rtol,
    atol=atol,
)


print("")
print("Initial Twiss functions:")
alpha_x, beta_x, alpha_y, beta_y, dispersion_x, dispersion_px = get_twiss(initial_beam)
print(f"  alpha_x={alpha_x:e} beta_x={beta_x:e} alpha_y={alpha_y:e} beta_y={beta_y:e}")
print(f"  dispersion_x={dispersion_x:e} dispersion_px={dispersion_px:e}")

atol = 0.0  # ignored
rtol = 3.5 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [alpha_x, beta_x, alpha_y, beta_y, dispersion_x],
    [
        1.340770e00,
        1.440253e01,
        -1.347747e00,
        4.602637e00,
        -2.667038e-01,
    ],
    rtol=rtol,
    atol=atol,
)

print("")
print("Final Twiss functions:")
alpha_x, beta_x, alpha_y, beta_y, dispersion_x, dispersion_px = get_twiss(final_beam)
print(f"  alpha_x={alpha_x:e} beta_x={beta_x:e} alpha_y={alpha_y:e} beta_y={beta_y:e}")
print(f"  dispersion_x={dispersion_x:e} dispersion_px={dispersion_px:e}")

atol = 0.0  # ignored
rtol = 3.5 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [beta_x, alpha_y, beta_y],
    [
        5.176642e00,
        -4.973010e00,
        4.687750e01,
    ],
    rtol=rtol,
    atol=atol,
)

# We use absolute tolerance for the following quantities, because these

assert np.allclose(
    [alpha_x],
    [
        0.0,
    ],
    atol=0.1,
)
assert np.allclose(
    [dispersion_x],
    [
        0.0,
    ],
    atol=3.0e-5,
)
