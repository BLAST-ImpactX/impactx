#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""The spin-orbit coupling of an RFCavity follows the Thomas-BMT equation.

The reference integrates the Lorentz force and the Thomas-BMT equation, with z as the
independent variable, through the paraxial field of the on-axis profile e(z):

    E_z = E0 e(z) cos(k t + phi),   E_x = -x/2 dE_z/dz,   B_y = x/2 dE_z/dt,

and likewise in y, with B_x = -y/2 dE_z/dt. The fields are normalized by q/(m c^2)
and q/(m c), in 1/m, and t is the time multiplied by c.

The transverse momenta of the RFCavity map are canonical momenta in the gauge without
a scalar potential: P_x = u_x + x/2 (E0/k) e'(z) sin(k t + phi), with the kinetic
momentum u_x = p_x/(m c).
"""

import numpy as np
import pytest

from impactx import Config, RefPart, elements

C_LIGHT = 299792458.0

if Config.precision == "SINGLE":
    rtol = 2.0e-3
else:
    rtol = 2.0e-5

CAVITIES = {
    "proton-2GeV": dict(
        species="proton",
        kin_energy_MeV=938.27208816,
        ds=0.7,
        escale=0.05,
        freq=7.0e8,
        phase=30.0,
        cos_coefficients=[0.4, 0.8, -0.2],
        sin_coefficients=[0.0, 0.3, 0.1],
    ),
    "proton-100MeV": dict(
        species="proton",
        kin_energy_MeV=100.0,
        ds=0.5,
        escale=0.02,
        freq=3.5e8,
        phase=-40.0,
        cos_coefficients=[1.0, 0.5, 0.2],
        sin_coefficients=[0.0, -0.2, 0.05],
    ),
    "electron-10MeV": dict(
        species="electron",
        kin_energy_MeV=10.0,
        ds=1.31879807,
        escale=8.0,
        freq=1.3e9,
        phase=60.0,
        cos_coefficients=[2.0, 0.3, -0.1],
        sin_coefficients=[0.0, 0.1, -0.03],
    ),
}


def make_cavity(cavity, mapsteps):
    ref = RefPart()
    ref.set_species(cavity["species"]).set_kin_energy_MeV(cavity["kin_energy_MeV"])
    el = elements.RFCavity(
        ds=cavity["ds"],
        escale=cavity["escale"],
        freq=cavity["freq"],
        phase=cavity["phase"],
        cos_coefficients=cavity["cos_coefficients"],
        sin_coefficients=cavity["sin_coefficients"],
        mapsteps=mapsteps,
    )
    return el, ref


def profile(z, cavity):
    """On-axis field profile e(z) and its slope e'(z) on 0 <= z <= ds"""
    ds = cavity["ds"]
    cos_coef = cavity["cos_coefficients"]
    sin_coef = cavity["sin_coefficients"]
    zc = z - 0.5 * ds
    e = 0.5 * cos_coef[0]
    ep = 0.0
    for j in range(1, len(cos_coef)):
        a = 2.0 * np.pi * j / ds
        e += cos_coef[j] * np.cos(a * zc) + sin_coef[j] * np.sin(a * zc)
        ep += a * (sin_coef[j] * np.cos(a * zc) - cos_coef[j] * np.sin(a * zc))
    return e, ep


def lorentz_tbmt(z, Y, cavity, G):
    """Derivative in z of the particle states Y: x, u_x, y, u_y, t, gamma and the
    spin rotation matrix"""
    E0 = cavity["escale"]
    k = 2.0 * np.pi * cavity["freq"] / C_LIGHT
    phi = np.deg2rad(cavity["phase"])

    x, ux, y, uy, t, gamma = Y[:, :6].T
    uz = np.sqrt(gamma**2 - 1.0 - ux**2 - uy**2)
    e, ep = profile(z, cavity)
    cos_term = np.cos(k * t + phi)
    sin_term = np.sin(k * t + phi)
    Ez = E0 * e * cos_term
    Ex = -0.5 * x * E0 * ep * cos_term
    Ey = -0.5 * y * E0 * ep * cos_term
    dEz_dt = -k * E0 * e * sin_term
    Bx = -0.5 * y * dEz_dt
    By = 0.5 * x * dEz_dt

    dY = np.empty_like(Y)
    dY[:, 0] = ux / uz
    dY[:, 1] = (gamma * Ex - uz * By) / uz
    dY[:, 2] = uy / uz
    dY[:, 3] = (gamma * Ey + uz * Bx) / uz
    dY[:, 4] = gamma / uz
    dY[:, 5] = (ux * Ex + uy * Ey + uz * Ez) / uz

    # Thomas-BMT: dS/dz = Omega x S
    u = np.stack([ux, uy, uz], axis=1)
    B = np.stack([Bx, By, np.zeros_like(x)], axis=1)
    E = np.stack([Ex, Ey, Ez], axis=1)
    u_dot_B = np.sum(u * B, axis=1)
    Omega = (
        -(
            (1.0 + G * gamma)[:, None] * B
            - (G / (gamma + 1.0) * u_dot_B)[:, None] * u
            - (G + 1.0 / (gamma + 1.0))[:, None] * np.cross(u, E)
        )
        / uz[:, None]
    )
    Ox, Oy, Oz = Omega.T
    zero = np.zeros_like(x)
    W = np.stack(
        [
            np.stack([zero, -Oz, Oy], axis=1),
            np.stack([Oz, zero, -Ox], axis=1),
            np.stack([-Oy, Ox, zero], axis=1),
        ],
        axis=1,
    )
    dY[:, 6:] = (W @ Y[:, 6:].reshape(-1, 3, 3)).reshape(-1, 9)
    return dY


def reference_maps(cavity, ref, nsteps=1000, delta=1.0e-6):
    """Linear map R and spin-orbit coupling A from central differences of
    particles integrated with RK4"""
    ds = cavity["ds"]
    E0 = cavity["escale"]
    k = 2.0 * np.pi * cavity["freq"] / C_LIGHT
    phi = np.deg2rad(cavity["phase"])
    G = ref.gyromagnetic_anomaly
    gamma0 = ref.gamma
    bg0 = ref.beta_gamma

    def gauge(z, Y):
        """P - u per unit transverse position"""
        _, ep = profile(z, cavity)
        return 0.5 * E0 / k * ep * np.sin(k * Y[:, 4] + phi)

    # the reference particle, then each coordinate shifted by +delta and -delta
    v = np.zeros((13, 6))
    for j in range(6):
        v[1 + 2 * j, j] = delta
        v[2 + 2 * j, j] = -delta
    Y = np.zeros((13, 15))
    Y[:, 0] = v[:, 0]
    Y[:, 2] = v[:, 2]
    Y[:, 4] = ref.t + v[:, 4]
    Y[:, 5] = gamma0 - bg0 * v[:, 5]
    Y[:, 1] = bg0 * v[:, 1] - gauge(0.0, Y) * Y[:, 0]
    Y[:, 3] = bg0 * v[:, 3] - gauge(0.0, Y) * Y[:, 2]
    Y[:, 6:] = np.eye(3).reshape(9)

    h = ds / nsteps
    for n in range(nsteps):
        z = n * h
        k1 = lorentz_tbmt(z, Y, cavity, G)
        k2 = lorentz_tbmt(z + 0.5 * h, Y + 0.5 * h * k1, cavity, G)
        k3 = lorentz_tbmt(z + 0.5 * h, Y + 0.5 * h * k2, cavity, G)
        k4 = lorentz_tbmt(z + h, Y + h * k3, cavity, G)
        Y = Y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    gamma1 = Y[0, 5]
    bg1 = np.sqrt(gamma1**2 - 1.0)
    w = np.empty((13, 6))
    w[:, 0] = Y[:, 0]
    w[:, 1] = (Y[:, 1] + gauge(ds, Y) * Y[:, 0]) / bg1
    w[:, 2] = Y[:, 2]
    w[:, 3] = (Y[:, 3] + gauge(ds, Y) * Y[:, 2]) / bg1
    w[:, 4] = Y[:, 4] - Y[0, 4]
    w[:, 5] = -(Y[:, 5] - gamma1) / bg1

    # the rotation vector of a small rotation
    S = Y[:, 6:].reshape(-1, 3, 3)
    lam = 0.5 * np.stack(
        [S[:, 2, 1] - S[:, 1, 2], S[:, 0, 2] - S[:, 2, 0], S[:, 1, 0] - S[:, 0, 1]],
        axis=1,
    )

    R = (w[1::2] - w[2::2]).T / (2.0 * delta)
    A = (lam[1::2] - lam[2::2]).T / (2.0 * delta)
    return R, A


def to_array(m, rows, cols):
    return np.array([[m[i, j] for j in range(1, cols + 1)] for i in range(1, rows + 1)])


@pytest.mark.parametrize("cavity", CAVITIES.values(), ids=CAVITIES.keys())
def test_spin_coupling_thomas_bmt(cavity):
    """Spin-orbit coupling and linear map against the Lorentz force and Thomas-BMT"""
    el, ref = make_cavity(cavity, mapsteps=2560)
    R = to_array(el.transfer_map(ref), 6, 6)
    A = to_array(el.spin_coupling, 3, 6)

    R_ref, A_ref = reference_maps(cavity, ref)

    # the reference model is the one of the map
    np.testing.assert_allclose(R, R_ref, rtol=0, atol=rtol * np.abs(R_ref).max())
    np.testing.assert_allclose(A, A_ref, rtol=0, atol=rtol * np.abs(A_ref).max())


def test_no_field_no_spin_rotation():
    """A cavity without a field does not couple the spin to the orbit"""
    cavity = dict(
        CAVITIES["proton-2GeV"], cos_coefficients=[0.0], sin_coefficients=[0.0]
    )
    el, ref = make_cavity(cavity, mapsteps=100)
    el.transfer_map(ref)
    np.testing.assert_array_equal(to_array(el.spin_coupling, 3, 6), 0.0)
