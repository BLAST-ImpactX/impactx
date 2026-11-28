#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Nilanjan Banerjee, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, elements

# a long proton bunch: Gauss2p5D assumes sigma_perp < beta*gamma*sigma_t
NUM_PARTICLES = 1000
SIGMA_X, SIGMA_Y, SIGMA_T = 1.0e-3, 2.0e-3, 1.0e-2  # rms bunch size in m
SHIFT = (2.0e-3, -1.0e-3, 5.0e-3)  # bunch offset from the reference particle in m


def _kicks(model, shift):
    """Momenta after one space-charge kick on a bunch at rest

    Every call creates the same macroparticles, translated as a whole by shift = (x, y, t).
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = model
    sim.diagnostics = False
    sim.slice_step_diagnostics = False

    try:
        sim.init_grids()

        ref = sim.beam.ref
        ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(250.0)
        qm_eev = ref.charge_qe / (ref.mass_MeV * 1.0e6)

        rng = np.random.default_rng(seed=42)
        x = rng.normal(0.0, SIGMA_X, NUM_PARTICLES) + shift[0]
        y = rng.normal(0.0, SIGMA_Y, NUM_PARTICLES) + shift[1]
        t = rng.normal(0.0, SIGMA_T, NUM_PARTICLES) + shift[2]
        zeros = np.zeros_like(x)
        sim.beam.add_n_particles(
            x, y, t, zeros, zeros, zeros, qm_eev, bunch_charge=1.0e-9
        )

        # one slice: a single space-charge kick, between two drift halves that keep the momenta
        sim.lattice.append(elements.Drift(name="drift", ds=1.0, nslice=1))
        sim.track_particles()

        # the same particle order in every run
        df = sim.beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
        momenta = ["momentum_x", "momentum_y", "momentum_t"]
        return {p: df[p].to_numpy().astype(np.float64) for p in momenta}
    finally:
        sim.finalize()


@pytest.mark.parametrize("model", ["Gauss3D", "Gauss2p5D"])
def test_gauss_space_charge_follows_the_bunch(model):
    """The analytic Gaussian space-charge fields are centered on the bunch centroid

    A bunch that is offset from the reference particle receives the same kicks as the
    same bunch on the reference orbit: there is no coherent self-kick.
    """
    centered = _kicks(model, shift=(0.0, 0.0, 0.0))
    shifted = _kicks(model, shift=SHIFT)

    # relative to the rms kick; in single precision, the Gauss2p5D line density
    # (binned and interpolated in t) amplifies the rounding of the shifted positions
    rtol = 1.0e-3 if Config.precision == "SINGLE" else 1.0e-10
    for p, kick in centered.items():
        rms = np.sqrt(np.mean(kick**2))
        assert rms > 0.0, f"{model}: no {p} kick"
        atol = rtol * rms
        if Config.precision == "SINGLE" and model == "Gauss3D" and p == "momentum_t":
            # FIXME: Gauss3D kicks at fixed t, and the transformation back to fixed s
            # rounds pt to multiples of ulp(gamma) / beta_gamma, here 1.5e-7 or 1.7e-3
            # of the rms kick. A particle of the shifted bunch can round to the
            # neighboring multiple.
            #   https://github.com/BLAST-ImpactX/impactx/issues/1705
            atol = 3.0e-3 * rms
        np.testing.assert_allclose(
            shifted[p], kick, rtol=0.0, atol=atol, err_msg=f"{model}: {p}"
        )
