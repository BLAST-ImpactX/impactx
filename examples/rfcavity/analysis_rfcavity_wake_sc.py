#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#

import glob
import re
from pathlib import Path

import numpy as np
import openpmd_api as io
import pandas as pd
import scipy.constants as const


# Load data from simulation
def read_time_series(file_pattern):
    """Read in all CSV files from each MPI rank (and potentially OpenMP
    thread). Concatenate into one Pandas dataframe.

    Returns
    -------
    pandas.DataFrame
    """

    def read_file(file_pattern):
        for filename in glob.glob(file_pattern):
            df = pd.read_csv(filename, delimiter=r"\s+")
            if "step" not in df.columns:
                step = int(re.findall(r"[0-9]+", filename)[0])
                df["step"] = step
            yield df

    return pd.concat(
        read_file(file_pattern),
        axis=0,
        ignore_index=True,
    )  # .set_index('id')


def data_is_double(file_pattern):
    """Detect float precision (single vs double) from a text diagnostic's digits."""
    text = "".join(Path(f).read_text() for f in glob.glob(file_pattern))
    return any(
        len(m.replace(".", "") if e else m.replace(".", "").strip("0")) >= 12
        for m, e in re.findall(r"(\d*\.\d+)([eE][+-]?\d+)?", text)
    )


# initial/final beam
series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
last_step = list(series.iterations)[-1]
initial = series.iterations[1].particles["beam"].to_df()
final = series.iterations[last_step].particles["beam"].to_df()

rbc_ref = read_time_series("diags/ref_particle.*")
is_double_ref = data_is_double("diags/ref_particle.*")

rbc = read_time_series("diags/reduced_beam_characteristics.*")
is_double = data_is_double("diags/reduced_beam_characteristics.*")

s = rbc["s"]
sigma_x = rbc["sigma_x"]
sigma_y = rbc["sigma_y"]
sigma_t = rbc["sigma_t"]
mean_pt = rbc["mean_pt"]
emittance_x = rbc["emittance_x"]
emittance_y = rbc["emittance_y"]
emittance_t = rbc["emittance_t"]

length = len(s) - 1

sig_xf = sigma_x.iloc[length]
sig_yf = sigma_y.iloc[length]
sig_tf = sigma_t.iloc[length]
mean_ptf = mean_pt.iloc[length]
emittance_xf = emittance_x.iloc[length]
emittance_yf = emittance_y.iloc[length]
emittance_tf = emittance_t.iloc[length]

s = rbc_ref["s"]
beta = rbc_ref["beta"]
gamma = rbc_ref["gamma"]

sf = s.iloc[-1]
betaf = beta.iloc[-1]
gammaf = gamma.iloc[-1]

# compare number of particles
num_particles = 1000000
assert num_particles == len(initial)
assert num_particles == len(final)

# IMPACT-Z: final beam moments
sigx_impactz = 0.356738e-03  # m
sigy_impactz = 0.246065e-03  # m
sigPhase_impactz = 0.109871e00  # deg
sigpx_impactz = 0.603754e-05  # rad
sigpy_impactz = 0.936262e-05  # rad
sigE_impactz = 0.623585e-01  # MeV
emitx_impactz = 0.100148e-05  # m-rad
emity_impactz = 0.100108e-05  # m-rad
emitPhaseE_impactz = 0.664925e-02  # deg-MeV
meanE_impactz = 0.279653e-01  # MeV (beam particles)

s_impactz = 0.593917e01
beta_impactz = 0.999999e00
gamma_impactz = 0.585160e03
ErestMeV = 0.51099895069
scale_frequency = 1.3e9
scale_wavelength = const.c / (2.0 * const.pi * scale_frequency)
p0_MeVc = beta_impactz * gamma_impactz * ErestMeV

sigt_predicted = sigPhase_impactz * scale_wavelength * const.pi / 180.0
sigpt_predicted = sigE_impactz / p0_MeVc
emitt_predicted = emitPhaseE_impactz * scale_wavelength * const.pi / (180.0 * ErestMeV)
meanpt_predicted = meanE_impactz / p0_MeVc

print("")
print("Reference particle path length and beam energy (final):")
print(f"  s={sf:e} gamma={gammaf:e}")

atol = 0.0  # ignored
rtol = 1.0e-5
if not is_double:
    # FIXME: single precision requires somewhat relaxed tolerance
    rtol = 1.0e-4
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sf, gammaf],
    [s_impactz, gamma_impactz],
    rtol=rtol,
    atol=atol,
)

print("")
print("Final Beam, predicted:")
print(f"  sigx={sigx_impactz:e} sigy={sigy_impactz:e} sigt={sigt_predicted:e}")
print(
    f"  emittance_xn={emitx_impactz:e} emittance_yn={emity_impactz:e} emittance_tn={emitt_predicted:e}"
)
print(f"  meanpt={meanpt_predicted:e}")
print("")
print("Final Beam:")
emittance_xn = emittance_xf * betaf * gammaf
emittance_yn = emittance_yf * betaf * gammaf
emittance_tn = emittance_tf * betaf * gammaf
print(f"  sigx={sig_xf:e} sigy={sig_yf:e} sigt={sig_tf:e}")
print(
    f"  emittance_xn={emittance_xn:e} emittance_yn={emittance_yn:e} emittance_tn={emittance_tn:e}"
)
print(f"  meanpt={mean_ptf:e}")

atol = 0.0  # ignored
rtol = 2.0 * num_particles**-0.5  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [sig_xf, sig_yf, sig_tf, emittance_xn, emittance_yn, emittance_tn],
    [
        sigx_impactz,
        sigy_impactz,
        sigt_predicted,
        emitx_impactz,
        emity_impactz,
        emitt_predicted,
    ],
    rtol=rtol,
    atol=atol,
)

atol = 0.0  # ignored
rtol = 0.02  # from random sampling of a smooth distribution
print(f"  rtol={rtol} (ignored: atol~={atol})")

assert np.allclose(
    [mean_ptf],
    [meanpt_predicted],
    rtol=rtol,
    atol=atol,
)
