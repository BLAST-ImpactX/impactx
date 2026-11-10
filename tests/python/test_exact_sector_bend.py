#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""Exact map of the sector bend, and its linear map, from strong to very weak bends.

In a sector bend of radius rc, the exact map moves a particle by a small amount
compared to rc. The elements ExactSbend and ExactCFbend (with only a dipole) must
reproduce this map to rounding, relative to the change of each coordinate, also for
|x/rc| << 1. Their linear maps must be accurate to rounding for any curvature and
gradient, including zero.
"""

from decimal import Decimal, localcontext

import numpy as np
import pytest

from impactx import Config, ImpactX, RefPart, elements

DS = 0.7  # element length [m]

# radius of curvature of the reference orbit [m]: a strong bend, a ring dipole,
# a reverse-bend quadrupole and a bend of 1e-12 degrees
RC = {
    "strong": 2.0,
    "ring": 40.7,
    "reverse": -1060.0,
    "1e-12deg": DS / np.deg2rad(1.0e-12),
}

# x [m], y [m], t [m], px, py, pt
PARTICLES = np.array(
    [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # reference orbit
        [0.0, -1.7e-3, 2.1e-3, 0.0, 0.0, 0.0],  # co-moving, offset in y and t
        [1.3e-3, 0.0, 0.0, 0.0, 0.0, 0.0],  # offset from the reference orbit
        [0.0, 0.0, 0.0, 1.0e-4, 0.0, 0.0],  # inclined
        [0.0, 0.0, 0.0, 0.0, 0.0, -2.0e-5],  # off-energy
        [1.3e-3, -1.7e-3, 2.1e-3, 1.1e-3, -0.8e-3, 1.0e-5],
        [-2.0e-3, 0.0, 0.0, -2.0e-3, 1.5e-3, -2.0e-5],
    ]
)

# dipoles of curvature 1/rc, which apply the exact map of the sector bend
ELEMENTS = {
    "ExactSbend": lambda rc: elements.ExactSbend(ds=DS, phi=np.rad2deg(DS / rc)),
    "ExactCFbend": lambda rc: elements.ExactCFbend(
        ds=DS, k_normal=[1.0 / rc], k_skew=[0.0]
    ),
}

COLUMNS = [
    "position_x",
    "position_y",
    "position_t",
    "momentum_x",
    "momentum_y",
    "momentum_t",
]

PREC = 60  # significant digits of the reference


def sin_cos(a):
    """sin(a) and cos(a) of a Decimal |a| <= pi, by their Taylor series"""
    tiny = Decimal(10) ** -(PREC + 5)
    s, c = Decimal(0), Decimal(0)
    term, k = Decimal(1), 0  # a^k / k!
    while k < 2 or abs(term) > tiny:
        if k % 4 == 0:
            c += term
        elif k % 4 == 1:
            s += term
        elif k % 4 == 2:
            c -= term
        else:
            s -= term
        k += 1
        term = term * a / k
    return s, c


def asin(v):
    """arcsin(v) of a Decimal |v| < 1, by the Taylor series of arctan"""
    # tan(asin(v)/4), with |z| <= tan(pi/8)
    z = v / (1 + (1 - v * v).sqrt())
    z = z / (1 + (1 + z * z).sqrt())
    tiny = Decimal(10) ** -(PREC + 5)
    total, power, k = Decimal(0), z, 0  # z^(2k+1)
    while abs(power) > tiny:
        total += (-1) ** k * power / (2 * k + 1)
        power *= z * z
        k += 1
    return 4 * total


def exact_sector_bend(beta, rc, ds, particle):
    """Exact map of the sector bend, evaluated with 60 significant digits

    D. L. Bruhwiler et al., Proc. EPAC 98, pp. 1171-1173 (1998)
    """
    with localcontext() as ctx:
        ctx.prec = PREC
        x, y, t, px, py, pt = (Decimal(float(v)) for v in particle)
        ibeta = 1 / Decimal(beta)
        rc = Decimal(rc)
        phi = Decimal(ds) / rc
        sin_phi, cos_phi = sin_cos(phi)

        pperp = (1 - 2 * pt * ibeta + pt * pt - py * py).sqrt()
        pzi = (pperp * pperp - px * px).sqrt()
        rho = rc + x
        pxf = px * cos_phi + (pzi - rho / rc) * sin_phi
        pzf = (pperp * pperp - pxf * pxf).sqrt()
        theta = phi + asin(px / pperp) - asin(pxf / pperp)

        out = (
            -rc + rho * cos_phi + rc * (pzf + px * sin_phi - pzi * cos_phi),
            y + theta * rc * py,
            t - theta * rc * (pt - ibeta) - phi * rc * ibeta,
            pxf,
            py,
            pt,
        )
        return [float(v) for v in out]


def coordinates(beam):
    """Particle coordinates in the order of insertion"""
    df = beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
    return df[COLUMNS].to_numpy(dtype=np.float64, copy=True)


def push(element, kin_energy_MeV, particles):
    """Push protons through an element

    Returns the stored initial coordinates, rounded to the precision of the build,
    the final coordinates and beta of the reference particle.
    """
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

    beam = sim.beam
    n = particles.shape[0]
    beam.add_n_particles(
        x=particles[:, 0],
        y=particles[:, 1],
        t=particles[:, 2],
        px=particles[:, 3],
        py=particles[:, 4],
        pt=particles[:, 5],
        qm=ref.qm_ratio_SI,
        w=np.ones(n),
    )

    initial = coordinates(beam)
    beta = ref.beta
    element.push(beam)
    final = coordinates(beam)
    sim.finalize()

    assert final.shape == (n, 6)
    return initial, final, beta


def check_exact_sector_bend(initial, final, beta, rc, ds):
    """Compare the final coordinates with the exact map of the sector bend"""
    expected = np.array([exact_sector_bend(beta, rc, ds, p) for p in initial])

    # py and pt do not change, and co-moving particles keep all coordinates
    np.testing.assert_array_equal(final[:, 4:6], initial[:, 4:6])
    comoving = np.all(initial[:, [0, 3, 4, 5]] == 0.0, axis=1)
    np.testing.assert_array_equal(final[comoving], initial[comoving])

    # the other particles: each coordinate to a few rounding errors of its size and of
    # its change, which is of order ds times the momentum deviation p_dev, including the
    # deflection h x of an offset particle; the arrival time t changes by ds p_dev / beta
    eps = np.finfo(np.float32 if Config.precision == "SINGLE" else np.float64).eps
    moving = ~comoving
    x0, y0, t0 = initial[moving, 0], initial[moving, 1], initial[moving, 2]
    px0, py0, pt0 = initial[moving, 3], initial[moving, 4], initial[moving, 5]
    p_dev = np.abs(px0) + np.abs(py0) + np.abs(pt0) / beta + np.abs(x0 / rc)
    scale = np.column_stack(
        [
            np.abs(x0) + ds * p_dev,
            np.abs(y0) + ds * p_dev,
            np.abs(t0) + ds * p_dev / beta,
            p_dev,
        ]
    )
    error = np.abs(final[moving, 0:4] - expected[moving, 0:4])
    assert np.all(error <= 64.0 * eps * scale), error / (eps * scale)


@pytest.mark.parametrize("element_name", ELEMENTS.keys())
@pytest.mark.parametrize("bend", RC.keys())
@pytest.mark.parametrize("kin_energy_MeV", [1.0e-3, 1.0e3], ids=["1keV", "1GeV"])
def test_exact_sector_bend(element_name, bend, kin_energy_MeV):
    """Dipoles reproduce the exact map of the sector bend, for any curvature."""
    rc = RC[bend]
    initial, final, beta = push(ELEMENTS[element_name](rc), kin_energy_MeV, PARTICLES)
    check_exact_sector_bend(initial, final, beta, rc, DS)


def test_exact_sector_bend_half_turn():
    """A half turn in one step, in which steep particles turn by more than a right angle
    relative to the reference orbit."""
    # x [m], y [m], t [m], px, py, pt
    particles = np.array(
        [
            [0.0, 0.0, 0.0, 0.75, 0.0, 0.0],
            [0.0, 0.0, 0.0, -0.75, 0.0, 0.0],
            [0.0, 0.0, 0.0, np.sqrt(0.5), 0.0, 0.0],  # turns by a right angle
            [1.3e-3, -1.7e-3, 2.1e-3, 0.72, -0.1, 1.0e-3],
            [0.0, 0.0, 0.0, 0.3, 0.0, 0.0],
        ]
    )
    rc, ds = 1.0, np.pi
    element = elements.ExactSbend(ds=ds, phi=180.0)
    initial, final, beta = push(element, 1.0e3, particles)
    check_exact_sector_bend(initial, final, beta, rc, ds)


# combined-function bends (curvature h [1/m], gradient k [1/m^2]): a weak dipole,
# focusing and defocusing in x, without curvature, and gx = h^2 + k = 0
CF_BENDS = {
    "dipole": (-1.0 / 1060.0, 0.0),
    "focusing": (1.0 / 40.7, 0.8),
    "defocusing": (-1.0 / 1060.0, -2.5),
    "straight": (0.0, 0.8),
    "gx=0": (0.5, -0.25),
}


def oscillator(g, s):
    """cos(w s), sin(w s)/w, (1 - cos(w s))/g and (sin(w s)/w - s)/g for g = w^2,
    and their hyperbolic counterparts for g < 0"""
    if g == 0:
        return Decimal(1), s, s * s / 2, -(s**3) / 6
    w = abs(g).sqrt()
    if g > 0:
        sin_ws, cos_ws = sin_cos(w * s)
    else:
        e = (w * s).exp()
        sin_ws, cos_ws = (e - 1 / e) / 2, (e + 1 / e) / 2
    return cos_ws, sin_ws / w, (1 - cos_ws) / g, (sin_ws / w - s) / g


def linear_sector_bend(beta, beta_gamma, h, k, ds):
    """Linear map of the sector bend with a gradient, evaluated with 60 digits

    Returns the map and the scale of the rounding error of R56: the size of its terms,
    which can cancel, where (beta gamma)^2 = gamma^2 - 1 of the reference particle's
    gamma is rounded relative to gamma^2, much larger at low energy.
    """
    with localcontext() as ctx:
        ctx.prec = PREC
        beta, beta_gamma = Decimal(beta), Decimal(beta_gamma)
        h, k, ds = Decimal(h), Decimal(k), Decimal(ds)
        gx, gy = k + h * h, -k
        cx, sx, vx, ux = oscillator(gx, ds)
        cy, sy, _, _ = oscillator(gy, ds)
        hb = h / beta

        R = [[Decimal(int(i == j)) for j in range(6)] for i in range(6)]
        R[0][0], R[0][1], R[0][5] = cx, sx, -vx * hb
        R[1][0], R[1][1], R[1][5] = -gx * sx, cx, -sx * hb
        R[2][2], R[2][3] = cy, sy
        R[3][2], R[3][3] = -gy * sy, cy
        R[4][0], R[4][1] = sx * hb, vx * hb
        R56_terms = (ds / beta_gamma**2, ux * hb * hb)
        R[4][5] = sum(R56_terms)
        gamma2_over_bg2 = (1 + beta_gamma**2) / beta_gamma**2
        return (
            np.array([[float(v) for v in row] for row in R]),
            float(abs(R56_terms[0]) * gamma2_over_bg2 + abs(R56_terms[1])),
        )


def check_linear_map(R, R_ref, R56_scale):
    """Each entry to a few rounding errors, R56 relative to the scale of its error"""
    eps = np.finfo(np.float32 if Config.precision == "SINGLE" else np.float64).eps
    scale = np.abs(R_ref)
    scale[4, 5] = R56_scale
    error = np.abs(R - R_ref)
    assert np.all(error <= 32.0 * eps * scale), error / (
        eps * np.maximum(scale, 1e-300)
    )


@pytest.mark.parametrize("bend", RC.keys())
@pytest.mark.parametrize("kin_energy_MeV", [1.0e-3, 1.0e3], ids=["1keV", "1GeV"])
def test_sbend_linear_map(bend, kin_energy_MeV):
    """The linear map of ExactSbend, for any curvature"""
    ref = RefPart()
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

    rc = RC[bend]
    R = ELEMENTS["ExactSbend"](rc).transfer_map(ref).to_numpy()
    R_ref, R56_scale = linear_sector_bend(ref.beta, ref.beta_gamma, 1.0 / rc, 0.0, DS)
    check_linear_map(R, R_ref, R56_scale)


@pytest.mark.parametrize("bend", CF_BENDS.keys())
@pytest.mark.parametrize("kin_energy_MeV", [1.0e-3, 1.0e3], ids=["1keV", "1GeV"])
def test_cfbend_linear_map(bend, kin_energy_MeV):
    """The linear map of ExactCFbend, for any curvature and gradient"""
    ref = RefPart()
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

    h, k = CF_BENDS[bend]
    element = elements.ExactCFbend(ds=DS, k_normal=[h, k], k_skew=[0.0, 0.0])
    R = element.transfer_map(ref).to_numpy()
    R_ref, R56_scale = linear_sector_bend(ref.beta, ref.beta_gamma, h, k, DS)
    check_linear_map(R, R_ref, R56_scale)
