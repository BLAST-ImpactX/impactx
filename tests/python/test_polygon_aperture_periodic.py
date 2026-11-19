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

import amrex.space3d as amr
from impactx import ImpactX, elements

# the rectangle x in [-0.3, 0.1] mm, y in [-0.1, 0.3] mm: not symmetric about x = 0 or y = 0
RECTANGLE_MM = ([-0.3, 0.1, 0.1, -0.3, -0.3], [-0.1, -0.1, 0.3, 0.3, -0.1])
# the same, widened to touch the edge of its period cell at x = -0.5 mm
TOUCHING_MM = ([-0.5, 0.1, 0.1, -0.5, -0.5], [-0.1, -0.1, 0.3, 0.3, -0.1])


def track_survivors(x, y, vertices_mm, **repeat):
    """Track particles at (x, y) through a repeated polygon.

    Returns the sorted (x, y) of the particles that pass on the MPI I/O rank, None on
    all other ranks.
    """
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False

    try:
        sim.init_grids()

        beam = sim.beam
        beam.ref.set_species("proton").set_kin_energy_MeV(1.0e3)
        # particles are added on the calling rank: add each one once
        io_rank = amr.ParallelDescriptor.IOProcessor()
        if io_rank:
            zero = np.zeros_like(x)
            qm_eev = beam.ref.charge_qe / (beam.ref.mass_MeV * 1.0e6)  # e / eV
            beam.add_n_particles(x, y, zero, zero, zero, zero, qm_eev, 1.0e-12)

        vertices_x, vertices_y = (np.array(c) * 1.0e-3 for c in vertices_mm)
        sim.lattice.append(
            elements.PolygonAperture(
                vertices_x=vertices_x, vertices_y=vertices_y, **repeat
            )
        )
        sim.track_particles()

        # collective: gathers the particles of all ranks on the I/O rank
        df = beam.to_df(
            local=False, root_rank=amr.ParallelDescriptor.IOProcessorNumber()
        )
        if not io_rank:
            return None
        if df is None or df.empty:
            return np.empty((0, 2))
        return np.array(sorted(df[["position_x", "position_y"]].to_numpy().tolist()))
    finally:
        sim.finalize()


@pytest.mark.parametrize(
    "vertices_mm, repeat, points_mm, expected_mm",
    [
        pytest.param(
            RECTANGLE_MM,
            dict(repeat_x=1.0e-3, repeat_y=1.0e-3),
            # inside the rectangle or a copy one period away, then two mirror images
            [(-0.2, 0.2), (-1.2, 0.2), (-0.2, -0.8), (0.2, 0.2), (-0.2, -0.2)],
            [(-1.2, 0.2), (-0.2, -0.8), (-0.2, 0.2)],
            id="repeat_xy",
        ),
        pytest.param(
            RECTANGLE_MM,
            dict(repeat_x=1.0e-3, repeat_y=1.0e-3, shift_odd_x=True),
            # copies in odd rows are shifted by half a period: inside a shifted copy in
            # rows 1 and -1 and inside a copy in row -2, then the unshifted positions
            # in rows 1 and -1 and a mirror image in row -2
            [
                (-0.6, 1.1),
                (-0.6, -0.9),
                (1.4, -0.9),
                (-0.2, -1.8),
                (-0.2, 1.1),
                (-0.2, -0.9),
                (0.2, -1.8),
            ],
            [(-0.6, -0.9), (-0.6, 1.1), (-0.2, -1.8), (1.4, -0.9)],
            id="shift_odd_x",
        ),
        pytest.param(
            TOUCHING_MM,
            dict(repeat_x=1.0e-3, repeat_y=1.0e-3),
            # inside the polygon at the touching edge, inside the copy beyond that edge
            # and inside another copy, then two points beside all copies
            [(-0.45, 0.2), (0.55, 0.2), (-1.3, -0.95), (0.3, 0.2), (-0.2, 0.45)],
            [(-0.45, 0.2), (0.55, 0.2), (-1.3, -0.95)],
            id="touching_cell_edge",
        ),
    ],
)
def test_periodic_polygon_aperture_is_not_mirrored(
    vertices_mm, repeat, points_mm, expected_mm
):
    """A periodic aperture repeats an asymmetric polygon by translation, not reflection."""
    x, y = (np.array(c) * 1.0e-3 for c in zip(*points_mm))
    survivors = track_survivors(x, y, vertices_mm, **repeat)

    if survivors is not None:
        expected = np.array(sorted(expected_mm)) * 1.0e-3
        assert survivors.shape == expected.shape
        np.testing.assert_allclose(survivors, expected, rtol=1.0e-6)
