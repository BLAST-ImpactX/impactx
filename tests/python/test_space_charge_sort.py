#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, ImpactX, elements


def track(space_charge, sort_interval):
    """Track the same beam with space charge, return its particles in memory order."""
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = space_charge
    sim.space_charge_sort_interval = sort_interval
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.max_level = 0
    sim.prob_relative = [1.1]
    if space_charge == "3D":
        sim.poisson_solver = "fft" if Config.have_fft else "multigrid"
        sim.n_cell = [16, 16, 16]
    else:
        sim.poisson_solver = "fft"
        sim.n_cell = [32, 32, 1]
        sim.blocking_factor_x = [16]
        sim.blocking_factor_y = [16]
        sim.blocking_factor_z = [1]
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(10.0)
    qm_eev = -1.0 / 0.510998950e6

    # the same particles in each run, independent of the random number generator of the simulation
    rng = np.random.default_rng(42)
    x, y, t = rng.normal(0.0, 1.0e-4, (3, 4096))
    zeros = np.zeros_like(x)
    sim.beam.add_n_particles(x, y, t, zeros, zeros, zeros, qm_eev, 1.0e-9)

    sim.lattice.append(elements.Drift(name="drift", ds=1.0, nslice=12))
    sim.track_particles()

    df = sim.beam.to_df(local=True).reset_index(drop=True)
    sim.finalize()
    return df


@pytest.mark.parametrize(
    "space_charge",
    [
        "3D",
        pytest.param(
            "2D", marks=pytest.mark.skipif(not Config.have_fft, reason="needs FFT")
        ),
    ],
)
def test_sorting_keeps_physics(space_charge):
    """Sorting changes the order of particles in memory, but not how they are tracked."""
    unsorted = track(space_charge, sort_interval=0)
    sorted_ = track(space_charge, sort_interval=5)

    # particle ids continue to count up between simulations: compare their order
    assert not np.array_equal(
        unsorted["idcpu"].to_numpy().argsort(), sorted_["idcpu"].to_numpy().argsort()
    )

    unsorted = unsorted.sort_values("idcpu")
    sorted_ = sorted_.sort_values("idcpu")

    # in single precision, the sum of the deposited charge depends more on the particle order
    tolerance = 1.0e-3 if unsorted["position_x"].dtype == np.float32 else 1.0e-9
    for attribute in unsorted.columns.drop("idcpu"):
        reference = unsorted[attribute].to_numpy()
        difference = reference - sorted_[attribute].to_numpy()
        assert np.abs(difference).max() <= tolerance * reference.std(), attribute


if __name__ == "__main__":
    test_sorting_keeps_physics("3D")
