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

import amrex.space3d as amr
from impactx import ImpactX, elements

# staggered slots: x in [-0.1, 0.1] mm around each hole, each row open over its height;
# inside slots in rows 0, 1, -1 and -2, then the unshifted position in row 1 and the
# shifted one in row 0
BRICK_POINTS = [(0.05, 0.45), (0.55, 0.6), (-0.45, -0.6), (1.03, -1.9)]
BRICK_POINTS += [(0.05, 0.6), (0.55, 0.45)]
BRICK_KEPT = [(0.05, 0.45), (0.55, 0.6), (-0.45, -0.6), (1.03, -1.9)]


def track_survivors(x, y, **aperture):
    """Track particles at (x, y) through a repeated Aperture.

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

        sim.lattice.append(elements.Aperture(**aperture))
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
    "aperture_mm, points_mm, expected_mm",
    [
        pytest.param(
            dict(aperture_x=0.15, aperture_y=0.10, shift_odd_x=True),
            # holes in rows 0, 1, -1 and 3 (odd rows shifted by half a period), then
            # the unshifted position in row 1 and two points beside the central hole
            [(0.0, 0.0), (0.5, 1.0), (-0.5, -1.05), (1.45, 2.95)]
            + [(0.0, 1.0), (0.3, 0.0), (0.0, 0.3)],
            [(0.0, 0.0), (0.5, 1.0), (-0.5, -1.05), (1.45, 2.95)],
            id="holes_shift_odd_x",
        ),
        pytest.param(
            dict(aperture_x=0.0, aperture_y=0.10, shift_odd_x=True),
            # no horizontal boundary: horizontal strips, the shift changes nothing
            [(0.37, 0.05), (1.81, 1.02), (-2.2, -0.95), (0.37, 0.3), (-1.4, 1.5)],
            [(0.37, 0.05), (1.81, 1.02), (-2.2, -0.95)],
            id="strips_shift_odd_x",
        ),
        pytest.param(
            dict(aperture_x=0.15, aperture_y=0.0),
            # no vertical boundary: vertical slits
            [(0.1, 0.7), (1.05, -2.3), (-0.9, 0.45), (0.5, 0.0), (-1.6, 1.2)],
            [(0.1, 0.7), (1.05, -2.3), (-0.9, 0.45)],
            id="slits",
        ),
        pytest.param(
            dict(
                aperture_x=0.30,
                aperture_y=0.20,
                shift_odd_x=True,
                shape="elliptical",
                action="absorb",
            ),
            # absorbing ellipses in rows 0, 1 and -1, then points outside all of them: the
            # unshifted position in row 1, inside the rectangle around an ellipse, and
            # between the rows
            [(0.0, 0.0), (0.5, 1.0), (-0.5, -1.0)]
            + [(0.0, 1.0), (0.25, 0.15), (0.5, 0.5)],
            [(0.0, 1.0), (0.25, 0.15), (0.5, 0.5)],
            id="ellipses_absorb_shift_odd_x",
        ),
        pytest.param(
            dict(aperture_x=0.10, aperture_y=0.50, shift_odd_x=True),
            # rows that touch
            BRICK_POINTS,
            BRICK_KEPT,
            id="touching_rows_shift_odd_x",
        ),
        pytest.param(
            dict(aperture_x=0.10, aperture_y=0.0, shift_odd_x=True),
            # no vertical boundary: each row is open over its height, the same slots
            BRICK_POINTS,
            BRICK_KEPT,
            id="open_rows_shift_odd_x",
        ),
        pytest.param(
            dict(aperture_x=0.10, aperture_y=0.70, shift_odd_x=True),
            # overlapping shifted rows: holes in rows 0 and 3, the copies of rows 1 and
            # -1 reaching into row 0, then two points beside all copies
            [(0.0, 0.0), (0.5, 0.4), (0.5, -0.4), (-0.5, 1.45)]
            + [(0.25, 0.0), (0.3, 0.45)],
            [(0.0, 0.0), (0.5, 0.4), (0.5, -0.4), (-0.5, 1.45)],
            id="overlapping_rows_shift_odd_x",
        ),
        pytest.param(
            dict(aperture_x=0.70, aperture_y=0.30, shape="elliptical"),
            # ellipses that overlap within their row: inside one or two of them, then
            # outside of all
            [(0.0, 0.0), (0.5, 0.0), (0.5, 1.1), (-1.45, -2.05)]
            + [(0.5, 0.25), (0.0, 0.5)],
            [(0.0, 0.0), (0.5, 0.0), (0.5, 1.1), (-1.45, -2.05)],
            id="overlapping_ellipses",
        ),
    ],
)
def test_periodic_aperture(aperture_mm, points_mm, expected_mm):
    """A repeated aperture acts as the union of its copies, one every 1 mm.

    A plane without a boundary fills its period.
    """
    aperture = {
        key: value * 1.0e-3 if key.startswith("aperture") else value
        for key, value in aperture_mm.items()
    }
    x, y = (np.array(c) * 1.0e-3 for c in zip(*points_mm))
    survivors = track_survivors(x, y, repeat_x=1.0e-3, repeat_y=1.0e-3, **aperture)

    if survivors is not None:
        expected = np.array(sorted(expected_mm)) * 1.0e-3
        assert survivors.shape == expected.shape
        np.testing.assert_allclose(survivors, expected, rtol=1.0e-6, atol=1.0e-12)
