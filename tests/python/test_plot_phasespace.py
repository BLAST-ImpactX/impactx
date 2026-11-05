#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""The values that ``beam.plot_phasespace()`` shows, read back from the figure."""

import importlib.util
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pytest

from impactx import ImpactX, amr

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("pandas") is None, reason="pandas is not available"
)

# half-widths of the beam in x, px, y, py, t, pt
#   powers of two stay exact in single precision, also in the units of the plot
HALF_WIDTHS = 2.0 ** np.array([-10, -12, -11, -13, -9, -8])

# the vertices of a box in phase space: all six coordinates are uncorrelated
SIGNS = np.array(list(itertools.product([-1.0, 1.0], repeat=6))).T

# weights, which differ between the vertices
WEIGHTS = 1.0e5 * np.prod(
    np.array([2.0, 5.0, 3.0, 7.0, 11.0, 13.0])[:, np.newaxis] ** (SIGNS > 0.0), axis=0
)


@pytest.fixture
def sim():
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.verbose = 0
    sim.init_grids()
    yield sim
    sim.finalize()


def add_beam(sim, species="electron", copies=1, **charge):
    """Add the beam, split into copies, on each MPI rank."""
    ref = sim.beam.ref
    ref.set_species(species).set_kin_energy_MeV(100.0)
    coordinates = np.repeat(SIGNS * HALF_WIDTHS[:, np.newaxis], copies, axis=1)
    x, px, y, py, t, pt = coordinates
    sim.beam.add_n_particles(
        x,
        y,
        t,
        px,
        py,
        pt,
        ref.qm_ratio_SI,
        **(charge or {"w": np.repeat(WEIGHTS / copies, copies)}),
    )


def panels(fig):
    """The axes of the x, y and t phase space plots, in this order."""
    panels = [ax for ax in fig.axes if ax.images]
    assert len(panels) == 3
    return panels


def bins(panel):
    """The values of the bins of a phase space plot, with the rows along the momentum."""
    return panel.images[0].get_array().filled(0.0)


def colorbar_label(panel):
    return panel.images[0].colorbar.ax.get_ylabel()


def expected_bins(plane, values):
    """Sum values over the 2x2 bins of a phase space plane, with rows along the momentum."""
    r, p = SIGNS[2 * plane], SIGNS[2 * plane + 1]
    return [
        [values[(r == r_sign) & (p == p_sign)].sum() for r_sign in (-1.0, 1.0)]
        for p_sign in (-1.0, 1.0)
    ]


@pytest.mark.parametrize("species", ["electron", "proton"])
@pytest.mark.parametrize("copies", [1, 2])
def test_histograms_show_charge(sim, species, copies):
    """The bins hold |q| w, however the beam is split into macroparticles."""
    add_beam(sim, species, copies)

    moments = sim.beam.beam_moments()
    fig = sim.beam.plot_phasespace(num_bins=2)
    if fig is None:  # the figure is on the MPI root rank
        return

    total_weight = WEIGHTS.sum() * amr.ParallelDescriptor.NProcs()
    charge = abs(sim.beam.ref.charge)
    assert abs(moments["charge_C"]) == pytest.approx(charge * total_weight, rel=1e-5)

    for plane, panel in enumerate(panels(fig)):
        np.testing.assert_allclose(
            bins(panel),
            np.multiply(
                expected_bins(plane, WEIGHTS),
                charge * amr.ParallelDescriptor.NProcs(),
            ),
            rtol=1e-5,
        )
        assert bins(panel).sum() == pytest.approx(abs(moments["charge_C"]), rel=1e-5)
        assert colorbar_label(panel) == "$|Q|$ [C/bin]"

    plt.close(fig)


def test_histograms_without_charge_show_macroparticles(sim):
    """Test particles of zero weight have no charge to show, but are still plotted."""
    add_beam(sim, bunch_charge=0.0)

    fig = sim.beam.plot_phasespace(num_bins=2)
    if fig is None:  # the figure is on the MPI root rank
        return

    for plane, panel in enumerate(panels(fig)):
        np.testing.assert_array_equal(
            bins(panel),
            np.multiply(
                expected_bins(plane, np.ones_like(WEIGHTS)),
                amr.ParallelDescriptor.NProcs(),
            ),
        )
        assert colorbar_label(panel) == "macroparticles/bin"

    plt.close(fig)


def test_axes_match_their_units(sim):
    """x, y, ct are in mm, px, py in mrad, and pt in units of p_0 c, without a factor."""
    add_beam(sim)

    fig = sim.beam.plot_phasespace(num_bins=2)
    if fig is None:  # the figure is on the MPI root rank
        return

    x, px, y, py, t, pt = HALF_WIDTHS
    for panel, r, p in zip(
        panels(fig), (x * 1e3, y * 1e3, t * 1e3), (px * 1e3, py * 1e3, pt)
    ):
        np.testing.assert_allclose(panel.images[0].get_extent(), [-r, r, -p, p])

    plt.close(fig)
