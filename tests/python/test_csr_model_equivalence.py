#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import numpy as np
import pytest
from scipy.constants import e as q_e

from impactx import Config, ImpactX, amr, elements, wakeconvolution

pytestmark = pytest.mark.skipif(
    not Config.have_fft, reason="The built-in CSR model requires ImpactX_FFT=ON"
)

# a reduced version of examples/chicane/run_chicane_csr.py
KIN_ENERGY_MEV = 5.0e3
BUNCH_CHARGE_C = 1.0e-9
NPART = 5000
CSR_BINS = 100
NSLICE = 10
QM_EEV = -1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV


def make_particles():
    """A deterministic, chirped electron bunch in s-coordinates.

    ImpactX distribution sampling is not reproducible across simulations in
    one process (the RNG state advances), so both runs load this same array.
    """
    rng = np.random.default_rng(42)
    x = rng.normal(0.0, 2.3e-5, NPART)
    y = rng.normal(0.0, 1.3e-5, NPART)
    t = rng.normal(0.0, 2.0e-4, NPART)
    px = rng.normal(0.0, 1.6e-6, NPART)
    py = rng.normal(0.0, 2.8e-6, NPART)
    # a strong longitudinal chirp, as in the chicane compression example
    pt = -36.0 * t + rng.normal(0.0, 1.0e-6, NPART)
    return x, y, t, px, py, pt


checked = {"convolution": False}


def analytic_csr(profile, context):
    """A Python CSR kick model replicating the built-in analytic model:
    dN/ds convolved with the steady-state CSR wake function."""
    lam = profile.charge.to_numpy(copy=True)
    n = profile.num_bins
    bin_size = profile.bin_size

    # dN/ds, the number density slope per bin (\see DerivativeCharge1D)
    slopes = np.diff(lam) / bin_size / q_e

    # steady-state CSR wake function on 2N support, in wrap-around order
    # (\see HandleWakefield.H)
    idx = np.arange(2 * n)
    s_wake = np.where(idx <= n, idx, idx - 2 * n) * bin_size
    wake = np.array([wakeconvolution.w_l_csr(s, context.rc, bin_size) for s in s_wake])
    wake[n] = 0.0

    # zero-padded FFT convolution, cropped to the first N bins
    # (\see convolve_fft)
    padded = np.zeros(2 * n)
    padded[:n] = slopes
    kick = np.fft.irfft(np.fft.rfft(padded) * np.fft.rfft(wake), n=2 * n)[:n] * bin_size

    if not checked["convolution"]:
        # cross-check the numpy convolution against the compiled one, once
        d_slopes = amr.PODVector_real_default.from_xp(slopes)
        d_wake = amr.PODVector_real_default.from_xp(wake)
        kick_ref = wakeconvolution.convolve_fft(d_slopes, d_wake, bin_size).to_numpy(
            copy=True
        )
        assert np.allclose(
            kick, kick_ref, rtol=1e-10, atol=1e-12 * np.max(np.abs(kick_ref))
        )
        checked["convolution"] = True

    return kick


def run_chicane(kick_model=None):
    """Track the reduced CSR chicane and return the final beam moments."""
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = False
    sim.csr = True
    sim.csr_bins = CSR_BINS
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    if kick_model is not None:
        sim.csr_kick_model = kick_model

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)

    x, y, t, px, py, pt = make_particles()
    sim.beam.add_n_particles(x, y, t, px, py, pt, QM_EEV, bunch_charge=BUNCH_CHARGE_C)

    ns = NSLICE
    rc = 10.3462283686195526  # bend radius (meters)
    psi = 0.048345620280243  # pole face rotation angle (radians)
    lb = 0.500194828041958  # bend arc length (meters)

    dr1 = elements.Drift(name="dr1", ds=5.0058489435, nslice=ns)
    dr2 = elements.Drift(name="dr2", ds=1.0, nslice=ns)
    dr3 = elements.Drift(name="dr3", ds=2.0, nslice=ns)
    sbend1 = elements.Sbend(name="sbend1", ds=lb, rc=-rc, nslice=ns)
    sbend2 = elements.Sbend(name="sbend2", ds=lb, rc=rc, nslice=ns)
    dipedge1 = elements.DipEdge(name="dipedge1", psi=-psi, rc=-rc, g=0.0, K2=0.0)
    dipedge2 = elements.DipEdge(name="dipedge2", psi=psi, rc=rc, g=0.0, K2=0.0)

    lattice_half = [sbend1, dipedge1, dr1, dipedge2, sbend2]
    sim.lattice.extend(lattice_half)
    sim.lattice.append(dr2)
    lattice_half.reverse()
    sim.lattice.extend(lattice_half)
    sim.lattice.append(dr3)

    sim.track_particles()

    moments = dict(sim.beam.beam_moments())
    sim.finalize()
    return moments


def test_csr_kick_model_matches_builtin():
    """An analytic Python kick model reproduces the built-in CSR model."""
    moments_builtin = run_chicane()
    moments_surrogate = run_chicane(kick_model=analytic_csr)

    assert checked["convolution"]

    rtol = 1.0e-8 if Config.precision == "DOUBLE" else 2.0e-4
    for key in [
        "sigma_x",
        "sigma_y",
        "sigma_t",
        "sigma_px",
        "sigma_py",
        "sigma_pt",
        "emittance_x",
        "emittance_y",
        "emittance_t",
    ]:
        assert math.isclose(
            moments_builtin[key], moments_surrogate[key], rel_tol=rtol
        ), (
            f"{key}: built-in {moments_builtin[key]} "
            f"vs surrogate {moments_surrogate[key]}"
        )
