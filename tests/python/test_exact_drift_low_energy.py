#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""Exact drift map of field-free exact elements, down to low energy.

The longitudinal momentum of the exact drift is
pz = sqrt(1 - 2 pt/beta + pt^2 - px^2 - py^2), so particles that move with the
reference particle (px = py = pt = 0) keep their coordinates, and the arrival time t
of other particles slips by ds * ((1/beta - pt)/pz - 1/beta).
"""

from decimal import Decimal, localcontext

import numpy as np
import pytest

from impactx import Config, ImpactX, elements

DS = 0.7  # element length [m]

# x [m], y [m], t [m], px, py, pt
PARTICLES = np.array(
    [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # reference orbit
        [1.3e-3, -1.7e-3, 2.1e-3, 0.0, 0.0, 0.0],  # co-moving, offset
        # from the origin, so that x, y, t at the exit are the displacements
        [0.0, 0.0, 0.0, 1.0e-4, 0.0, 0.0],  # slightly inclined: slippage ~ px^2
        [0.0, 0.0, 0.0, 1.1e-3, -0.8e-3, 1.0e-5],  # inclined, off-energy
        [0.0, 0.0, 0.0, -2.0e-3, 1.5e-3, -2.0e-5],
    ]
)

# elements in field-free configurations, which apply the exact drift map
ELEMENTS = {
    "ExactDrift": lambda: elements.ExactDrift(ds=DS),
    "ExactQuad": lambda: elements.ExactQuad(ds=DS, k=0.0),
    "ExactMultipole": lambda: elements.ExactMultipole(
        ds=DS, k_normal=[0.0, 0.0], k_skew=[0.0, 0.0]
    ),
    "ExactSbend": lambda: elements.ExactSbend(ds=DS, phi=0.0),
    "ExactCFbend": lambda: elements.ExactCFbend(ds=DS, k_normal=[0.0], k_skew=[0.0]),
}

COLUMNS = [
    "position_x",
    "position_y",
    "position_t",
    "momentum_x",
    "momentum_y",
    "momentum_t",
]


def exact_drift(beta, ds, particle):
    """Exact drift map, evaluated with 50 significant digits"""
    with localcontext() as ctx:
        ctx.prec = 50
        x, y, t, px, py, pt = (Decimal(float(v)) for v in particle)
        ibeta = 1 / Decimal(beta)
        ds = Decimal(ds)
        pz = (1 - 2 * pt * ibeta + pt * pt - px * px - py * py).sqrt()
        out = (
            x + ds * px / pz,
            y + ds * py / pz,
            t - ds * (ibeta + (pt - ibeta) / pz),
            px,
            py,
            pt,
        )
        return [float(v) for v in out]


def coordinates(beam):
    """Particle coordinates in the order of insertion"""
    df = beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
    return df[COLUMNS].to_numpy(dtype=np.float64, copy=True)


@pytest.mark.parametrize("element_name", ELEMENTS.keys())
@pytest.mark.parametrize("kin_energy_MeV", [1.0e-3, 1.0e3], ids=["1keV", "1GeV"])
def test_exact_drift_low_energy(element_name, kin_energy_MeV):
    """Field-free exact elements reproduce the exact drift map.

    At 1 keV, the proton beta is 1.46e-3, so 1/beta^2 and 1/(beta*gamma)^2
    are 4.7e5 and differ by exactly 1.
    """
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

    beam = sim.beam
    n = PARTICLES.shape[0]
    beam.add_n_particles(
        x=PARTICLES[:, 0],
        y=PARTICLES[:, 1],
        t=PARTICLES[:, 2],
        px=PARTICLES[:, 3],
        py=PARTICLES[:, 4],
        pt=PARTICLES[:, 5],
        qm=ref.qm_ratio_SI,
        w=np.ones(n),
    )

    # the stored values, rounded to the precision of the build
    initial = coordinates(beam)
    beta = ref.beta

    element = ELEMENTS[element_name]()
    ds = element.ds
    element.push(beam)

    final = coordinates(beam)
    sim.finalize()

    assert final.shape == (n, 6)
    expected = np.array([exact_drift(beta, ds, p) for p in initial])

    # the momenta do not change, and co-moving particles keep all coordinates
    np.testing.assert_array_equal(final[:, 3:6], initial[:, 3:6])
    comoving = np.all(initial[:, 3:6] == 0.0, axis=1)
    np.testing.assert_array_equal(final[comoving], initial[comoving])

    # displacements of the other particles: x and y to a few rounding errors of their
    # size, and the small slippage of t to a few rounding errors of the terms of dt/ds,
    # of order pt/beta^2, pt and p_perp^2/beta, which cancel at high energy
    eps = np.finfo(np.float32 if Config.precision == "SINGLE" else np.float64).eps
    moving = ~comoving
    np.testing.assert_allclose(
        final[moving, 0:2], expected[moving, 0:2], rtol=32.0 * eps, atol=0.0
    )
    px, py, pt = initial[moving, 3], initial[moving, 4], initial[moving, 5]
    dt_terms = ds * (np.abs(pt) * (1.0 + 1.0 / beta**2) + (px**2 + py**2) / beta)
    np.testing.assert_allclose(
        (final[moving, 2] - expected[moving, 2]) / dt_terms,
        0.0,
        rtol=0.0,
        atol=32.0 * eps,
    )
