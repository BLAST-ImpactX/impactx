#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
from analysis_APL_ChrPlasmaLens import get_beams, get_moments, get_twiss

# initial/final beam
(initial, beam_final, final) = get_beams()

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

(betax, betay, alphax, alphay) = get_twiss(initial)
print(f"  betax={betax}[m],betay={betay}[m],alphax={alphax},alphay={alphay}")

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

# Compare to analytical values
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
sigx, sigy, sigt, emittance_x, emittance_y, emittance_t = get_moments(final)
s_ref = beam_final.get_attribute("s_ref")
gamma_ref = beam_final.get_attribute("gamma_ref")
print(f"  sigx={sigx:e} sigy={sigy:e} sigt={sigt:e}")
print(
    f"  emittance_x={emittance_x:e} emittance_y={emittance_y:e} emittance_t={emittance_t:e}\n"
    f"  s_ref={s_ref:e} gamma_ref={gamma_ref:e}"
)

(betax, betay, alphax, alphay) = get_twiss(final)
print(f"  betax={betax}[m],betay={betay}[m],alphax={alphax},alphay={alphay}")

atol = 0.0  # ignored
rtol = 2.2 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

# Compare to analytical values
assert np.allclose(
    [sigx, sigy, sigt, emittance_x, emittance_y, emittance_t, s_ref, gamma_ref],
    [
        2.737665020201518e-05,
        2.737665020201518e-05,
        0.001,
        2.548491664266332e-08,
        2.548491664266332e-08,
        1e-06,
        20e-3,
        3.923902e02,
    ],
    rtol=rtol,
    atol=atol,
)
