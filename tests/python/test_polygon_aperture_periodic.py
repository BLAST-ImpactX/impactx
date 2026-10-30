#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Nikita Kuklev
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, elements


def test_periodic_polygon_aperture_is_not_mirrored():
    """A periodic aperture repeats an asymmetric polygon by translation, not reflection."""
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    beam = sim.beam
    beam.ref.set_species("proton").set_kin_energy_MeV(1.0e3)
    # (x, y) in mm: inside the polygon or a copy one period away, then two mirror images
    x = np.array([-0.2, -1.2, -0.2, 0.2, -0.2]) * 1.0e-3
    y = np.array([0.2, 0.2, -0.8, 0.2, -0.2]) * 1.0e-3
    zero = np.zeros_like(x)
    beam.add_n_particles(x, y, zero, zero, zero, zero, beam.ref.qm_ratio_SI, 1.0e-12)

    # the rectangle x in [-0.3, 0.1] mm, y in [-0.1, 0.3] mm, repeated every 1 mm
    sim.lattice.append(
        elements.PolygonAperture(
            vertices_x=[-0.3e-3, 0.1e-3, 0.1e-3, -0.3e-3, -0.3e-3],
            vertices_y=[-0.1e-3, -0.1e-3, 0.3e-3, 0.3e-3, -0.1e-3],
            repeat_x=1.0e-3,
            repeat_y=1.0e-3,
        )
    )
    sim.track_particles()
    df = beam.to_df(local=True)
    survivors = (
        [] if df is None else df[["position_x", "position_y"]].to_numpy().tolist()
    )

    sim.finalize()

    expected = np.array([[-1.2, 0.2], [-0.2, -0.8], [-0.2, 0.2]]) * 1.0e-3
    assert len(survivors) == len(expected)
    np.testing.assert_allclose(sorted(survivors), expected, rtol=1.0e-6)
