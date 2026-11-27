#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Ji Qiang, Axel Huebl
# License: BSD-3-Clause-LBNL
#
"""The 2D expanding beam of analysis_expanding_fft_2D.py without longitudinal extent.

The 2D space-charge model of an unbunched beam does not depend on the longitudinal
coordinate, so a beam whose particles all have the same time of flight must expand
in the same way: the size doubles over the drift and the emittance stays zero.
"""

import numpy as np
import openpmd_api as io
from scipy.stats import moment

series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

num_particles = 10000
assert num_particles == len(initial)
assert num_particles == len(final)

# tolerances from random sampling of a smooth distribution, as in analysis_expanding_fft_2D.py
for name, beam, sig_expected, rtol in (
    ("Initial", initial, 5.0e-4, 1.5 * num_particles**-0.5),
    ("Final", final, 1.0e-3, 3.5 * num_particles**-0.5),
):
    beam = beam.select_dtypes(np.floating).astype(np.float64)
    sigx = moment(beam["position_x"], moment=2) ** 0.5
    sigy = moment(beam["position_y"], moment=2) ** 0.5
    sigt = moment(beam["position_t"], moment=2) ** 0.5
    cov = beam.cov(ddof=0)
    emittance_x = (
        cov["position_x"]["position_x"] * cov["momentum_x"]["momentum_x"]
        - cov["position_x"]["momentum_x"] ** 2
    ) ** 0.5
    emittance_y = (
        cov["position_y"]["position_y"] * cov["momentum_y"]["momentum_y"]
        - cov["position_y"]["momentum_y"] ** 2
    ) ** 0.5
    print(
        f"{name} Beam: sigx={sigx:e} sigy={sigy:e} sigt={sigt:e} "
        f"emittance_x={emittance_x:e} emittance_y={emittance_y:e} rtol={rtol:e}"
    )
    assert np.allclose([sigx, sigy], [sig_expected, sig_expected], rtol=rtol, atol=0.0)
    assert sigt == 0.0
    assert emittance_x < 1.0e-8
    assert emittance_y < 1.0e-8
