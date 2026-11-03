#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, distribution, twiss

npart = 100_000

# spins are unit vectors up to rounding
atol_norm = 1.0e-12 if Config.precision == "DOUBLE" else 1.0e-5
# sampling error of a mean spin component or second moment: at most
# 1/sqrt(3 npart) = 2e-3, which dominates the error of the inverse Langevin
# approximation that sets the concentration of the vMF distribution
atol_moments = 1.0e-2


def sample_spins(polarization):
    """Spins of a 1 GeV electron beam with the given polarization vector."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.spin = True
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(1.0e3)
    distr = distribution.Waterbag(
        **twiss(
            beta_x=1.0,
            beta_y=1.0,
            beta_t=1.0,
            emitt_x=1.0e-9,
            emitt_y=1.0e-9,
            emitt_t=1.0e-6,
        )
    )
    sim.add_particles(1.0e-9, distr, npart, distribution.SpinvMF(*polarization))

    df = sim.beam.to_df(local=True)
    spins = df[["spin_x", "spin_y", "spin_z"]].to_numpy(dtype=np.float64)

    sim.finalize()
    return spins


@pytest.mark.parametrize(
    "polarization",
    [(0.0, 0.0, 0.0), (0.0, 0.0, 0.5), (0.0, 0.0, -0.5), (0.3, -0.2, 0.4)],
    ids=["unpolarized", "along_z", "along_minus_z", "oblique"],
)
def test_spin_polarization(polarization):
    """Spins are unit vectors whose mean is the polarization vector."""
    spins = sample_spins(polarization)

    assert len(spins) == npart
    assert np.isfinite(spins).all()
    np.testing.assert_allclose(np.linalg.norm(spins, axis=1), 1.0, atol=atol_norm)
    np.testing.assert_allclose(spins.mean(axis=0), polarization, atol=atol_moments)


def test_spin_unpolarized_isotropic():
    """Without polarization, spins are uniform on the unit sphere:
    <s_i s_j> = delta_ij / 3."""
    spins = sample_spins((0.0, 0.0, 0.0))

    second_moments = spins.T @ spins / len(spins)
    np.testing.assert_allclose(second_moments, np.eye(3) / 3.0, atol=atol_moments)
