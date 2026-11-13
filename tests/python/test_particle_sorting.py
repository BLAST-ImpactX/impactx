#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Andrew Myers
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
Sorting the beam particles before space charge deposition (``algo.sort_intervals``)
must only change the order in which particles are stored, not the physics:
the beam moments after tracking must agree with an unsorted run up to round-off.
"""

import numpy as np
import pytest

import amrex.space3d as amr
from impactx import Config, ImpactX, elements

KIN_ENERGY_MEV = 250.0
BUNCH_CHARGE_C = 1.0e-9
NPART = 10000
QM_EEV = -1.0 / 0.510998950 / 1e6  # electron charge/mass in e / eV


def _deterministic_beam():
    """A fixed, reproducible beam, so every run starts from bit-identical particles."""
    rng = np.random.default_rng(seed=42)

    x = rng.normal(0.0, 4.472135955e-4, NPART)
    y = rng.normal(0.0, 4.472135955e-4, NPART)
    t = rng.normal(0.0, 9.12241869e-7, NPART)
    px = rng.normal(0.0, 1.0e-4, NPART)
    py = rng.normal(0.0, 1.0e-4, NPART)
    pt = np.zeros(NPART)

    return x, y, t, px, py, pt


def _run(sort_intervals, for_deposition=False, idx_type=(0, 0, 0), bin_size=(1, 1, 1)):
    """Track the beam through a drift with 3D space charge, four solves in total.

    All sorting options are set explicitly in every run, because the parameter
    database is shared by all simulations in this process.

    :return: the beam moments after tracking and the particle ids in storage order
    """
    sim = ImpactX()

    sim.n_cell = [16, 16, 16]
    sim.particle_shape = 2
    sim.space_charge = "3D"
    sim.dynamic_size = True
    sim.prob_relative = [3.0]
    sim.slice_step_diagnostics = False
    sim.diagnostics = False

    sim.sort_intervals = sort_intervals
    sim.sort_particles_for_deposition = for_deposition
    sim.sort_idx_type = list(idx_type)
    sim.sort_bin_size = list(bin_size)

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)
    ref.z = 0.0

    if amr.ParallelDescriptor.IOProcessor():
        x, y, t, px, py, pt = _deterministic_beam()
        sim.beam.add_n_particles(
            x, y, t, px, py, pt, QM_EEV, bunch_charge=BUNCH_CHARGE_C
        )

    sim.lattice.extend([elements.Drift(name="d1", ds=2.0, nslice=4)])

    try:
        sim.track_particles()

        # the options read back as set
        assert sim.sort_intervals == sort_intervals
        assert sim.sort_particles_for_deposition == for_deposition
        assert sim.sort_idx_type == list(idx_type)
        assert sim.sort_bin_size == list(bin_size)

        moments = sim.beam.beam_moments()
        df = sim.beam.to_df(local=True)
        ids = np.array([]) if df is None else df["idcpu"].to_numpy()
    finally:
        sim.finalize()

    return moments, ids


# sorting for deposition is only implemented in AMReX for these backends
HAVE_DEPOSITION_SORT = Config.gpu_backend in ("CUDA", "HIP")
needs_deposition_sort = pytest.mark.skipif(
    not HAVE_DEPOSITION_SORT, reason="sorting for deposition needs CUDA or HIP"
)


@pytest.mark.parametrize(
    "for_deposition,idx_type,bin_size",
    [
        pytest.param(True, (0, 0, 0), (1, 1, 1), marks=needs_deposition_sort),
        pytest.param(True, (2, 2, 2), (1, 1, 1), marks=needs_deposition_sort),
        (False, (0, 0, 0), (1, 1, 1)),
        (False, (0, 0, 0), (2, 2, 2)),
    ],
    ids=["deposition-cell", "deposition-mixed", "bin-1", "bin-2"],
)
def test_sorting_preserves_physics(for_deposition, idx_type, bin_size):
    """Sorted runs reorder particles but track them like the unsorted run."""
    ref_moments, ref_ids = _run(sort_intervals=-1)
    moments, ids = _run(
        sort_intervals=1,
        for_deposition=for_deposition,
        idx_type=idx_type,
        bin_size=bin_size,
    )

    # the same particles survive, but they are stored in a different order
    #   ids continue counting across simulations in one process, so compare the
    #   storage order by id rank within each run instead of the raw ids
    assert len(ids) == len(ref_ids)
    if amr.ParallelDescriptor.NProcs() == 1:
        ranks = np.argsort(np.argsort(ids))
        ref_ranks = np.argsort(np.argsort(ref_ids))
        assert not np.array_equal(ranks, ref_ranks), "particles were not sorted"

    # charge deposition sums in a different order: agreement up to round-off
    rtol = 1e-4 if Config.precision == "SINGLE" else 1e-10
    for key in ["sig_x", "sig_y", "sig_t", "sig_px", "sig_py", "sig_pt", "emittance_x"]:
        assert np.isclose(moments[key], ref_moments[key], rtol=rtol, atol=0.0), (
            f"{key}: sorted {moments[key]!r} != unsorted {ref_moments[key]!r}"
        )


def test_sort_invalid_bin_size():
    """Sorting in bins smaller than one cell is an error."""
    with pytest.raises(RuntimeError, match="sort_bin_size"):
        _run(sort_intervals=1, for_deposition=False, bin_size=(0, 1, 1))


@needs_deposition_sort
def test_sort_invalid_idx_type():
    """Sorting to a grid type other than 0, 1 or 2 is an error."""
    with pytest.raises(RuntimeError, match="sort_idx_type"):
        _run(sort_intervals=1, for_deposition=True, idx_type=(3, 0, 0))


@pytest.mark.skipif(HAVE_DEPOSITION_SORT, reason="backend sorts for deposition")
def test_sort_for_deposition_unavailable():
    """Without CUDA or HIP, sorting for deposition is a clear error, not an abort."""
    with pytest.raises(RuntimeError, match="only implemented for the CUDA and HIP"):
        _run(sort_intervals=1, for_deposition=True)
