#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

import amrex.space3d as amr
from impactx import ImpactX, elements


def expected_ez0_transport_map(ds, bz, beta_gamma):
    """Return the analytic ChrAcc transport map for zero Ez."""
    alpha = 0.5 * bz
    drift_factor = ds / beta_gamma

    matrix = np.eye(6)
    matrix[4, 5] = ds / beta_gamma**2

    if alpha == 0.0:
        matrix[0, 1] = ds
        matrix[2, 3] = ds
        return matrix

    theta = alpha * drift_factor
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    sin2 = sin_theta**2
    cos2 = cos_theta**2
    sin_cos = sin_theta * cos_theta

    matrix[0, 0] = cos2
    matrix[0, 1] = beta_gamma * sin_cos / alpha
    matrix[0, 2] = sin_cos
    matrix[0, 3] = beta_gamma * sin2 / alpha
    matrix[1, 0] = -alpha * sin_cos / beta_gamma
    matrix[1, 1] = cos2
    matrix[1, 2] = -alpha * sin2 / beta_gamma
    matrix[1, 3] = sin_cos
    matrix[2, 0] = -sin_cos
    matrix[2, 1] = -beta_gamma * sin2 / alpha
    matrix[2, 2] = cos2
    matrix[2, 3] = beta_gamma * sin_cos / alpha
    matrix[3, 0] = alpha * sin2 / beta_gamma
    matrix[3, 1] = -sin_cos
    matrix[3, 2] = -alpha * sin_cos / beta_gamma
    matrix[3, 3] = cos2

    return matrix


@pytest.mark.parametrize("bz_scale", [0.0, 1.0], ids=["no-bz", "bz"])
def test_element_push(bz_scale):
    """
    This tests using ImpactX without a lattice.
    """
    sim = ImpactX()

    # set numerical parameters and IO control
    sim.space_charge = False
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    # basic beam parameters
    kin_energy_MeV = 2000.0  # reference energy (kinetic)

    # set reference particle
    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
    qm_eev = ref.charge_qe / (ref.mass_MeV * 1.0e6)  # electron charge/mass in e / eV
    initial_ref = ref.copy()

    # set test particles
    pc = sim.beam
    npart = 20

    ptmin = -3.0e-3
    ptmax = 3.0e-3

    # Add off-axis particles in the field-free case so that the transverse
    # contribution to the zero-Ez t update is exercised. The finite-Bz case
    # remains on axis so that the inverse ChrDrift isolates the t update.
    if amr.ParallelDescriptor.IOProcessor():
        dpt = np.linspace(ptmin, ptmax, npart)
        zero_arr = np.zeros(npart)
        if bz_scale == 0.0:
            x = np.linspace(-1.0e-3, 1.0e-3, npart)
            y = np.linspace(0.75e-3, -0.75e-3, npart)
            px = np.linspace(-2.0e-4, 2.0e-4, npart)
            py = np.linspace(1.5e-4, -1.5e-4, npart)
        else:
            x = y = px = py = zero_arr
        pc.add_n_particles(
            x,
            y,
            zero_arr,
            px,
            py,
            dpt,
            qm_eev,
            bunch_charge=0.0,
        )

    # store initial particles
    df_initial = pc.to_df()

    # design the accelerator lattice
    ns = 1  # number of slices per ds in the element

    # add beam diagnostics
    monitor = elements.BeamMonitor("monitor", backend="h5")

    # bend radius (> 0)
    ks_value = 1.0

    # length for this test should be one period
    ds_value = 2.0 * np.pi / (ref.gyromagnetic_anomaly * ks_value)

    idrift_chr = elements.ChrDrift(name="idrift_chr", ds=-ds_value, nslice=ns)

    bz_value = bz_scale * ks_value * ref.beta_gamma
    ez_value = 0.0
    sol_chr = elements.ChrAcc(
        name="sol_chr", ds=ds_value, ez=ez_value, bz=bz_value, nslice=ns
    )
    np.testing.assert_allclose(
        sol_chr.transfer_map(initial_ref).to_numpy(),
        expected_ez0_transport_map(ds=ds_value, bz=bz_value, beta_gamma=ref.beta_gamma),
        atol=1.0e-12,
        rtol=0.0,
    )

    # set the lattice
    sim.lattice.append(monitor)
    sim.lattice.append(sol_chr)
    sim.lattice.append(idrift_chr)
    sim.lattice.append(monitor)

    # run simulation
    sim.track_particles()

    pc = sim.beam
    df_final = pc.to_df()

    PHASE_COLS = [
        "position_x",
        "position_y",
        "position_t",
        "momentum_x",
        "momentum_y",
        "momentum_t",
    ]

    REF_COLS = ["x", "y", "z", "t", "px", "py", "pz", "pt", "s"]

    for c in PHASE_COLS:
        np.testing.assert_allclose(
            df_final[c].to_numpy(),
            df_initial[c].to_numpy(),
            atol=1.0e-12,
            rtol=0,
            err_msg=f"Roundtrip mismatch in {c}",
        )
    for c in REF_COLS:
        np.testing.assert_allclose(
            getattr(sim.beam.ref, c),
            getattr(initial_ref, c),
            atol=1.0e-12,
            rtol=0,
            err_msg=f"Reference-particle roundtrip mismatch in {c}",
        )

    # clean shutdown
    sim.finalize()
