#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Nikita Kuklev
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, elements


@pytest.mark.parametrize("h, ks", [(0.3, 0.02), (-0.3, 0.02), (-0.3, 0.0)])
def test_exact_cfbend_skew_dipole(h, ks):
    """The skew dipole field ks/(1 + h x) kicks an on-axis particle by ks * ds."""
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    beam = sim.beam
    beam.ref.set_species("proton").set_kin_energy_MeV(1.0e3)
    beam.add_n_particles(
        [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], beam.ref.qm_ratio_SI, 1.0e-12
    )

    ds = 0.7
    elements.ExactCFbend(ds=ds, k_normal=[h], k_skew=[ks]).push(beam)
    df = beam.to_df(local=True)
    phase_space = df[
        [
            "position_x",
            "position_y",
            "position_t",
            "momentum_x",
            "momentum_y",
            "momentum_t",
        ]
    ].to_numpy()

    sim.finalize()

    rtol = 1.0e-6 if Config.precision == "SINGLE" else 1.0e-14
    assert np.isfinite(phase_space).all()
    assert phase_space[0, 4] == pytest.approx(ks * ds, rel=rtol, abs=0.0)
