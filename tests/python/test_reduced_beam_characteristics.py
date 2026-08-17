#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from pathlib import Path

import pytest

from impactx import ImpactX, distribution, elements


def test_reduced_beam_characteristics_momentum_columns():
    """
    The min/max momentum columns of the reduced beam characteristics hold the
    momentum extrema, not a repeat of the mean.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = True
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(
        1.0e-9,
        distribution.Waterbag(
            lambdaX=3.9984884770e-5,
            lambdaY=3.9984884770e-5,
            lambdaT=1.0e-3,
            lambdaPx=2.6623538760e-5,
            lambdaPy=2.6623538760e-5,
            lambdaPt=2.0e-3,
            muxpx=-0.846574929020762,
            muypy=0.846574929020762,
            mutpt=0.0,
        ),
        1000,
    )

    sim.lattice.append(elements.Drift(name="d1", ds=0.25))
    sim.track_particles()

    moments = sim.beam.beam_moments()
    sim.finalize()

    diag_file = next(Path("diags").glob("reduced_beam_characteristics_final.*"))
    header, *rows = diag_file.read_text().splitlines()
    columns = header.split()
    values = dict(zip(columns, (float(v) for v in rows[-1].split())))

    for coordinate in ("px", "py", "pt"):
        mean = values[f"mean_{coordinate}"]
        minimum = values[f"min_{coordinate}"]
        maximum = values[f"max_{coordinate}"]

        assert minimum < mean < maximum
        assert minimum == pytest.approx(moments[f"min_{coordinate}"])
        assert maximum == pytest.approx(moments[f"max_{coordinate}"])
