#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import re

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


def sample_radius(u, cdf, r):
    """The radius that Thermal samples for uniform draws u in [0, 1)."""
    off = np.maximum(np.searchsorted(cdf, u, side="left") - 1, 0)
    tv = (u - cdf[off]) / (cdf[off + 1] - cdf[off])
    return r[off] + tv * (r[off + 1] - r[off])


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
    "normalize_halo",
    [0.08195, 0.06],
    ids=["above_1", "below_1"],
)
def test_bithermal_radial_cdf(normalize_halo):
    """
    The radial CDFs rise monotonically from 0 to exactly 1, also if the radial
    densities integrate to more (0.08195) or less (0.06) than 1.
    """
    sim = thermal_sim()

    distr = distribution.Thermal(
        k=k,
        kT=kT,
        kT_halo=900.0e-6,
        normalize=0.54226,
        normalize_halo=normalize_halo,
        halo=0.05,
    )
    r, *cdfs = distr.radial_cdf(bunch_charge_C, sim.beam.ref)
    sim.finalize()

    r = np.asarray(r)
    assert r[0] > 0.0
    assert np.all(np.diff(r) > 0.0)

    for cdf in cdfs:
        cdf = np.asarray(cdf)
        nbins = len(cdf) - 1
        assert len(r) == len(cdf)
        assert cdf[0] == 0.0
        assert cdf[nbins] == 1.0
        assert np.all(np.diff(cdf) >= 0.0)

        # the radial lookup of the largest uniform draw below 1 stays in the table
        u = np.nextafter(1.0, 0.0)
        off = max(0, np.searchsorted(cdf, u, side="left") - 1)
        assert off + 1 <= nbins


def test_thermal_radial_cdf_charge_sign():
    """Electrons and positrons have the same radial density."""
    sim = thermal_sim()

    distr = distribution.Thermal(
        k=k,
        kT=kT,
        kT_halo=kT,
        normalize=0.41604661,
        normalize_halo=0.0,
        halo=0.0,
    )
    cdfs = {}
    for species in ["electron", "positron"]:
        sim.beam.ref.set_species(species).set_kin_energy_MeV(kin_energy_MeV)
        cdfs[species] = distr.radial_cdf(bunch_charge_C, sim.beam.ref)
    sim.finalize()

    core = np.asarray(cdfs["electron"][1])
    assert np.isfinite(core).all()
    assert core[-1] == 1.0
    assert cdfs["electron"] == cdfs["positron"]


def test_thermal_radial_cdf_radius():
    """The radii of the sampled particles follow the tabulated radial CDF."""
    sim = thermal_sim()

    distr = distribution.Thermal(
        k=k,
        kT=kT,
        kT_halo=kT,
        normalize=0.41604661,
        normalize_halo=0.0,
        halo=0.0,
    )
    r, core, _ = (np.asarray(v) for v in distr.radial_cdf(bunch_charge_C, sim.beam.ref))
    sim.add_particles(bunch_charge_C, distr, npart)

    # radius in the beam rest frame
    df = sim.beam.to_df(local=True)
    bg = sim.beam.ref.beta_gamma
    r_particles = np.sqrt(
        df["position_x"] ** 2 + df["position_y"] ** 2 + (bg * df["position_t"]) ** 2
    )
    sim.finalize()

    quantiles = np.array([0.25, 0.5, 0.75])
    assert np.allclose(
        np.quantile(r_particles, quantiles),
        sample_radius(quantiles, core, r),
        rtol=0.03,
    )


@pytest.mark.parametrize(
    "population, weight, normalize, normalize_halo",
    [("core", "0.95", 0.0, 0.08195), ("halo", "0.05", 0.54226, 0.0)],
    ids=["core", "halo"],
)
def test_bithermal_population_without_density(
    population, weight, normalize, normalize_halo
):
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
    message = (
        f"the {population} population is a fraction of {weight} of the beam, "
        "but its radial density integrates to 0."
    )
    with pytest.raises(RuntimeError, match=re.escape(message)):
        sim.add_particles(bunch_charge_C, distr, npart)
    assert sim.beam.total_number_of_particles() == 0

    sim.finalize()


@pytest.mark.parametrize(
    "normalize, normalize_halo, halo, parameters",
    [
        (0.5, 0.0, 0.0, "normalize and kT"),
        (0.54226, 0.1, 0.05, "normalize, normalize_halo, kT and kT_halo"),
    ],
    ids=["thermal", "bithermal"],
)
def test_diverging_radial_density(normalize, normalize_halo, halo, parameters):
    """
    A normalization above its self-consistent value makes the radial density
    diverge. Core and halo share one potential, so the error names the
    parameters of all populations in the beam.
    """
    sim = thermal_sim()

    distr = distribution.Thermal(
        k=k,
        kT=kT,
        kT_halo=900.0e-6,
        normalize=normalize,
        normalize_halo=normalize_halo,
        halo=halo,
    )
    message = f"the radial density diverges.*Please check {re.escape(parameters)}\\.$"
    with pytest.raises(RuntimeError, match=message):
        sim.add_particles(bunch_charge_C, distr, npart)
    assert sim.beam.total_number_of_particles() == 0

    sim.finalize()
