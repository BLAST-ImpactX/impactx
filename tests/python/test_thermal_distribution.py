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

from impactx import ImpactX, distribution

# thermal and bithermal beams of the EPAC2004 benchmarks,
# see examples/epac2004_benchmarks/
kin_energy_MeV = 0.1
bunch_charge_C = 1.4285714285714285714e-10
k = 6.283185307179586
kT = 36.0e-6
npart = 10_000


def thermal_sim():
    """A 0.1 MeV proton reference particle without any particles."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
    return sim


def test_thermal_without_halo():
    """
    A thermal beam has no halo population: its halo density is zero, but it is
    never sampled.
    """
    sim = thermal_sim()

    distr = distribution.Thermal(
        k=k,
        kT=kT,
        kT_halo=kT,
        normalize=0.41604661,
        normalize_halo=0.0,
        halo=0.0,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    df = sim.beam.to_df(local=True)
    assert len(df) == npart
    phase_space = df.filter(regex="^(position|momentum)_").to_numpy()
    assert phase_space.shape[1] == 6
    assert np.isfinite(phase_space).all()

    sim.finalize()


@pytest.mark.parametrize(
    "population, normalize, normalize_halo",
    [("core", 0.0, 0.08195), ("halo", 0.54226, 0.0)],
    ids=["core", "halo"],
)
def test_bithermal_population_without_density(population, normalize, normalize_halo):
    """
    A population that is part of the beam cannot be sampled without a radial
    density.
    """
    sim = thermal_sim()

    distr = distribution.Thermal(
        k=k,
        kT=kT,
        kT_halo=900.0e-6,
        normalize=normalize,
        normalize_halo=normalize_halo,
        halo=0.05,
    )
    with pytest.raises(RuntimeError, match=f"the {population} population"):
        sim.add_particles(bunch_charge_C, distr, npart)
    assert sim.beam.total_number_of_particles() == 0

    sim.finalize()
