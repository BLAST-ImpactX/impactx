#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
from analysis_APL import read_time_series

rbc = read_time_series("diags/reduced_beam_characteristics.*")

print("Initial Beam:")

sigx = rbc["sigma_x"].iloc[0]
sigy = rbc["sigma_y"].iloc[0]
sigt = rbc["sigma_t"].iloc[0]
emittance_x = rbc["emittance_x"].iloc[0]
emittance_y = rbc["emittance_y"].iloc[0]
emittance_t = rbc["emittance_t"].iloc[0]

print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}"
)

betax = rbc["beta_x"].iloc[0]
betay = rbc["beta_y"].iloc[0]
alphax = rbc["alpha_x"].iloc[0]
alphay = rbc["alpha_y"].iloc[0]

print(f"  betax={betax}[m],betay={betay}[m],alphax={alphax},alphay={alphay}")

atol = 0.0  # ignored
rtol = 1e-5
print(f"  rtol={rtol} (ignored: atol~={atol})")

# Compare initial beam to analytical values
assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t],
    [
        2.737665020201518e-05,
        2.737665020201518e-05,
        0.001,
        2.548491664266332e-08,
        2.548491664266332e-08,
        1e-06,
    ],
    rtol=rtol,
    atol=atol,
)

print("")
print("Final Beam:")

sigx = rbc["sigma_x"].iloc[-1]
sigy = rbc["sigma_y"].iloc[-1]
sigt = rbc["sigma_t"].iloc[-1]
emittance_x = rbc["emittance_x"].iloc[-1]
emittance_y = rbc["emittance_y"].iloc[-1]
emittance_t = rbc["emittance_t"].iloc[-1]

s_ref = rbc["s"].iloc[-1]

print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}\n"
    f"  s_ref={s_ref:e}"
)

betax = rbc["beta_x"].iloc[-1]
betay = rbc["beta_y"].iloc[-1]
alphax = rbc["alpha_x"].iloc[-1]
alphay = rbc["alpha_y"].iloc[-1]
print(f"  betax={betax}[m],betay={betay}[m],alphax={alphax},alphay={alphay}")

atol = 0.0  # ignored
rtol = 1e-5
print(f"  rtol={rtol} (ignored: atol~={atol})")

# Compare final beam to analytical values
assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, s_ref],
    [
        2.737665020201518e-05,
        2.737665020201518e-05,
        0.001,
        2.548491664266332e-08,
        2.548491664266332e-08,
        1e-06,
        20e-3,
    ],
    rtol=rtol,
    atol=atol,
)
