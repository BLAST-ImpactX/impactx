#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

import amrex.space3d as amr
from impactx import Config, ImpactX, elements

# magnetic shielding of the proton in H-
# T. Kilich and K. Pachucki, Phys. Rev. A 112, 042813 (2025)
SIGMA = 24.243189e-6

# the anomaly of H- is stored with 10 decimals
RTOL = 1.0e-10 if Config.precision != "SINGLE" else 1.0e-5


def spin_angle(species):
    """Spin rotation in rad of an on-axis particle in a longitudinal magnetic field."""
    c = 299792458.0  # speed of light in m/s
    field = 1.0  # longitudinal magnetic field in T
    ds = 1.0  # element length in m
    gamma = 2.0

    sim = ImpactX()

    sim.space_charge = False
    sim.spin = True
    sim.slice_step_diagnostics = False

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species(species)
    ref.set_kin_energy_MeV((gamma - 1.0) * ref.mass_MeV)

    if amr.ParallelDescriptor.IOProcessor():
        zero = np.zeros(1)
        one = np.ones(1)
        sim.beam.add_n_particles(
            zero,
            zero,
            zero,
            zero,
            zero,
            zero,
            ref.charge_qe / (ref.mass_MeV * 1.0e6),
            w=one,
            sx=one,
            sy=zero,
            sz=zero,
        )

    # the field is normalized by q / (m c)
    bz = ref.qm_ratio_SI * field / c
    sim.lattice.append(elements.ChrAcc(name="sol", ds=ds, ez=0.0, bz=bz))
    sim.track_particles()

    df = sim.beam.to_df()
    angle = float(np.arctan2(df["spin_y"].iloc[0], df["spin_x"].iloc[0]))

    # Thomas-BMT equation in a longitudinal magnetic field
    np.testing.assert_allclose(
        angle,
        -(1.0 + ref.gyromagnetic_anomaly) * bz * ds / ref.beta_gamma,
        rtol=RTOL,
        atol=0.0,
    )

    sim.finalize()

    return angle


def test_hminus_spin_precession():
    """
    The closed-shell electrons of H- carry no spin: its spin and magnetic moment are
    those of the shielded proton, irrespective of the charge of the ion.

    At equal gamma, a longitudinal magnetic field rotates the spin by an angle that is
    proportional to the magnetic moment. H- and the proton thus rotate in the same
    direction, with a ratio of (1 - sigma).
    """
    angle_proton = spin_angle("proton")
    angle_hminus = spin_angle("Hminus")

    np.testing.assert_allclose(
        angle_hminus,
        angle_proton * (1.0 - SIGMA),
        rtol=RTOL,
        atol=0.0,
    )
