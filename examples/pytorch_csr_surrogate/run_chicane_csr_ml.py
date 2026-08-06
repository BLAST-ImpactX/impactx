#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Auralee Edelen, Chris Mayes
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
#
# Train a small neural network surrogate for the steady-state CSR kick on the
# fly and couple it into the chicane example as a user-defined CSR kick model
# (sim.csr_kick_model), following the approach of
# A. L. Edelen et al., "Neural Network Solver for Coherent Synchrotron
# Radiation Wakefield Calculations in Accelerator-Based Charged Particle
# Beams", in Proc. IPAC'22, WEPOMS013 (https://arxiv.org/abs/2203.07542).
#
# The network maps the binned, normalized longitudinal charge profile to the
# dimensionless per-bin CSR kick shape. The physical kick follows from the
# exact scaling of the steady-state model,
#     kick = Q * kappa(R) / (q_e * bin_size^(4/3)) * G(shape),
# so one trained network generalizes over bunch charge, bunch length
# (compression stage), and bend radius.

import sys

import numpy as np
from scipy.constants import e as q_e
from scipy.constants import epsilon_0

from impactx import Config, ImpactX, distribution, elements

try:
    import torch
    from torch import nn
except ImportError:
    print("Warning: Cannot import PyTorch. Skipping test.")
    sys.exit(42)  # ImpactX special return code for skipped tests

if not Config.have_openpmd:
    print("Warning: ImpactX was built without openPMD support. Skipping test.")
    sys.exit(42)

# reproducible training
torch.manual_seed(7)
rng = np.random.default_rng(7)

CSR_BINS = 150  # matches sim.csr_bins below


def kappa(R):
    """CSR wake strength kappa = q_e^2 / (2 pi epsilon_0 3^(1/3) R^(2/3)) [J m^(1/3)].

    This is the prefactor of the bin-integrated steady-state CSR wake
    function, see w_l_csr in ImpactX (Saldin et al., NIMA 398, 373 (1997)).
    """
    return q_e**2 / (2.0 * np.pi * epsilon_0 * 3.0 ** (1.0 / 3.0) * R ** (2.0 / 3.0))


def normalized_csr_wake(n):
    """The dimensionless CSR wake on 2n support in wrap-around order.

    w_l_csr(s = j * bin_size, R, bin_size) = -kappa / bin_size^(1/3) * wn[j]
    with wn independent of R and bin_size.
    """
    j = np.arange(2 * n)
    s = np.where(j <= n, j, j - 2 * n).astype(np.float64)
    wn = 1.5 * (
        np.heaviside(s + 0.5, 1.0) * np.abs(s + 0.5) ** (2.0 / 3.0)
        - np.heaviside(s - 0.5, 1.0) * np.abs(s - 0.5) ** (2.0 / 3.0)
    )
    wn[n] = 0.0
    return wn


def analytic_g(p, wn):
    """Dimensionless per-bin kick shape G for a normalized profile p.

    p has n + 1 entries (per-bin charge fractions, sum ~ 1); the result has
    n entries. This replicates the built-in pipeline (charge derivative,
    zero-padded FFT convolution with the wake) in normalized units:
    kick = Q * kappa(R) / (q_e * bin_size^(4/3)) * G.
    """
    n = len(p) - 1
    padded = np.zeros(2 * n)
    padded[:n] = np.diff(p)
    return -np.fft.irfft(np.fft.rfft(padded) * np.fft.rfft(wn), n=2 * n)[:n]


def self_test_normalization():
    """Verify the normalized pipeline against the physical one for one profile."""
    n = CSR_BINS
    R = 10.3462283686195526
    sigma = 2.0e-5  # bunch length [m]
    charge = 1.0e-9  # bunch charge [C]

    z_min, z_max = -4.0 * sigma, 4.0 * sigma
    bin_size = (z_max - z_min) / (n - 1)
    z = z_min + np.arange(n + 1) * bin_size
    lam = np.exp(-0.5 * (z / sigma) ** 2)
    lam *= charge / (np.sum(lam) * bin_size)  # line density [C/m]

    # physical pipeline, as in the built-in model
    slopes = np.diff(lam) / bin_size / q_e
    j = np.arange(2 * n)
    s = np.where(j <= n, j, j - 2 * n).astype(np.float64) * bin_size
    wake = (
        -1.5
        * kappa(R)
        / bin_size
        * (
            np.heaviside(s + bin_size / 2, 1.0)
            * np.abs(s + bin_size / 2) ** (2.0 / 3.0)
            - np.heaviside(s - bin_size / 2, 1.0)
            * np.abs(s - bin_size / 2) ** (2.0 / 3.0)
        )
    )
    wake[n] = 0.0
    padded = np.zeros(2 * n)
    padded[:n] = slopes
    kick_physical = (
        np.fft.irfft(np.fft.rfft(padded) * np.fft.rfft(wake), n=2 * n)[:n] * bin_size
    )

    # normalized pipeline
    p = lam * bin_size / charge
    f_scale = charge * kappa(R) / (q_e * bin_size ** (4.0 / 3.0))
    kick_normalized = f_scale * analytic_g(p, normalized_csr_wake(n))

    assert np.allclose(
        kick_physical,
        kick_normalized,
        rtol=1.0e-10,
        atol=1.0e-12 * np.max(np.abs(kick_physical)),
    ), "normalized CSR pipeline does not match the physical one"


def sample_shapes(m, n):
    """Sample m normalized profiles (n + 1 bins) as 1-2 Gaussian mixtures."""
    shapes = np.zeros((m, n + 1))
    xi = np.arange(n + 1) / (n - 1)  # normalized bin coordinate
    for i in range(m):
        n_peaks = rng.integers(1, 3)
        weights = rng.uniform(0.3, 1.0, n_peaks)
        p = np.zeros(n + 1)
        for k in range(n_peaks):
            mu = rng.uniform(0.3, 0.7)
            sig = rng.uniform(0.08, 0.25)
            p += weights[k] * np.exp(-0.5 * ((xi - mu) / sig) ** 2)
        shapes[i] = p / np.sum(p)
    return shapes


def train_model():
    """Train a small MLP mapping profile shape -> dimensionless kick shape."""
    n = CSR_BINS
    wn = normalized_csr_wake(n)

    n_train, n_test = 800, 100
    shapes = sample_shapes(n_train + n_test, n)
    targets = np.array([analytic_g(p, wn) for p in shapes])
    g_scale = targets.std()

    X = torch.from_numpy(shapes)
    Y = torch.from_numpy(targets / g_scale)
    X_train, X_test = X[:n_train], X[n_train:]
    Y_train, Y_test = Y[:n_train], Y[n_train:]

    model = nn.Sequential(
        nn.Linear(n + 1, 64),
        nn.Tanh(),
        nn.Linear(64, 64),
        nn.Tanh(),
        nn.Linear(64, n),
    ).double()

    optimizer = torch.optim.Adam(model.parameters(), lr=2.0e-3)
    loss_fn = nn.MSELoss()
    for epoch in range(500):
        optimizer.zero_grad()
        loss = loss_fn(model(X_train), Y_train)
        loss.backward()
        optimizer.step()
        if epoch % 100 == 0:
            print(f"  epoch {epoch:4d}: loss={loss.item():.3e}")

    with torch.no_grad():
        pred = model(X_test)
        rel_err = (pred - Y_test).norm(dim=1) / Y_test.norm(dim=1)
    rel_err = rel_err.numpy()
    print(
        f"  test kick shape error: median={np.median(rel_err):.1%} "
        f"max={np.max(rel_err):.1%}"
    )
    assert np.median(rel_err) < 0.1, "NN CSR surrogate did not train well enough"

    return model, g_scale


def make_csr_kick_nn(model, g_scale):
    """Wrap the trained network as an ImpactX CSR kick model."""

    def csr_kick_nn(profile, context):
        lam = profile.charge.to_numpy(copy=True)
        bin_size = profile.bin_size

        charge = np.sum(lam) * bin_size  # bunch charge [C]
        p = lam * bin_size / charge  # normalized profile (sum ~ 1)
        f_scale = charge * kappa(context.rc) / (q_e * bin_size ** (4.0 / 3.0))

        with torch.no_grad():
            g = model(torch.from_numpy(p)).numpy()

        return f_scale * g_scale * g  # per-bin longitudinal force [N]

    return csr_kick_nn


if __name__ == "__main__":
    self_test_normalization()

    # multi-threaded training, before AMReX/OpenMP starts
    torch.set_num_threads(min(4, torch.get_num_threads()))
    print("Training the NN CSR surrogate ...")
    model, g_scale = train_model()

    # PyTorch's threaded defaults interfere with AMReX OpenMP during
    # tracking (see impactx#773, pyamrex#322)
    torch.set_num_threads(1)

    # the chicane example with the NN CSR kick model,
    # following examples/chicane/run_chicane_csr.py
    sim = ImpactX()

    sim.particle_shape = 2  # B-spline order
    sim.space_charge = False
    sim.csr = True
    sim.csr_bins = CSR_BINS
    sim.csr_kick_model = make_csr_kick_nn(model, g_scale)
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    # load a 5 GeV electron beam with an initial
    # normalized transverse rms emittance of 1 um
    kin_energy_MeV = 5.0e3  # reference energy
    bunch_charge_C = 1.0e-9  # used with space charge
    npart = 10000  # number of macro particles

    #   reference particle
    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

    #   particle bunch
    distr = distribution.Gaussian(
        lambdaX=2.2951017632e-5,
        lambdaY=1.3084093142e-5,
        lambdaT=5.5555553e-8,
        lambdaPx=1.598353425e-6,
        lambdaPy=2.803697378e-6,
        lambdaPt=2.000000000e-6,
        muxpx=0.933345606203060,
        muypy=0.933345606203060,
        mutpt=0.999999961419755,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    # add beam diagnostics
    monitor = elements.BeamMonitor("monitor", backend="h5")

    # design the accelerator lattice
    ns = 25  # number of slices per ds in the element
    rc = 10.3462283686195526  # bend radius (meters)
    psi = 0.048345620280243  # pole face rotation angle (radians)
    lb = 0.500194828041958  # bend arc length (meters)

    # Drift elements
    dr1 = elements.Drift(name="dr1", ds=5.0058489435, nslice=ns)
    dr2 = elements.Drift(name="dr2", ds=1.0, nslice=ns)
    dr3 = elements.Drift(name="dr3", ds=2.0, nslice=ns)

    # Bend elements
    sbend1 = elements.Sbend(name="sbend1", ds=lb, rc=-rc, nslice=ns)
    sbend2 = elements.Sbend(name="sbend2", ds=lb, rc=rc, nslice=ns)

    # Dipole Edge Focusing elements
    dipedge1 = elements.DipEdge(name="dipedge1", psi=-psi, rc=-rc, g=0.0, K2=0.0)
    dipedge2 = elements.DipEdge(name="dipedge2", psi=psi, rc=rc, g=0.0, K2=0.0)

    lattice_half = [sbend1, dipedge1, dr1, dipedge2, sbend2]
    # assign a segment with the first half of the lattice
    sim.lattice.append(monitor)
    sim.lattice.extend(lattice_half)
    sim.lattice.append(dr2)
    lattice_half.reverse()
    # extend the lattice by a reversed half
    sim.lattice.extend(lattice_half)
    sim.lattice.append(dr3)
    sim.lattice.append(monitor)

    # run simulation
    sim.track_particles()

    # clean shutdown
    sim.finalize()
