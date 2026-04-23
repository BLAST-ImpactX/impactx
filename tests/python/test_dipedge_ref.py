#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import numpy as np

import amrex.space3d as amr
from impactx import Config, ImpactX, elements


def _track_particle_dipedge(modify_ref_part, edge_angle, g, rc):
    sim = ImpactX()

    # set numerical parameters and IO control
    sim.space_charge = False
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    # load initial beam parameters
    kin_energy_MeV = 0.8e3  # reference energy
    bunch_charge_C = 1.0e-9  # used with space charge
    npart = 10000  # number of macro particles

    #   reference particle
    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)
    qm_eev = 1.0 / 938.27208816 / 1e6  # electron charge/mass in e / eV

    beam = sim.beam

    # define a single particle that initially coincides with reference orbit
    dx = [0.0]
    dpx = [0.0]
    dy = [0.0]
    dpy = [0.0]
    dt = [0.0]
    dpt = [0.0]
    if not Config.have_gpu:  # initialize using cpu-based PODVectors
        dx_podv = amr.PODVector_real_std()
        dy_podv = amr.PODVector_real_std()
        dt_podv = amr.PODVector_real_std()
        dpx_podv = amr.PODVector_real_std()
        dpy_podv = amr.PODVector_real_std()
        dpt_podv = amr.PODVector_real_std()
    else:  # initialize on device using arena/gpu-based PODVectors
        dx_podv = amr.PODVector_real_arena()
        dy_podv = amr.PODVector_real_arena()
        dt_podv = amr.PODVector_real_arena()
        dpx_podv = amr.PODVector_real_arena()
        dpy_podv = amr.PODVector_real_arena()
        dpt_podv = amr.PODVector_real_arena()

    for p_dx in dx:
        dx_podv.push_back(p_dx)
    for p_dy in dy:
        dy_podv.push_back(p_dy)
    for p_dt in dt:
        dt_podv.push_back(p_dt)
    for p_dpx in dpx:
        dpx_podv.push_back(p_dpx)
    for p_dpy in dpy:
        dpy_podv.push_back(p_dpy)
    for p_dpt in dpt:
        dpt_podv.push_back(p_dpt)

    beam.add_n_particles(
        dx_podv, dy_podv, dt_podv, dpx_podv, dpy_podv, dpt_podv, qm_eev, bunch_charge_C
    )

    # design the accelerator lattice)
    ns = 1  # number of slices per ds in the element

    dipedge1 = elements.DipEdge(
        name="dipedge1",
        psi=edge_angle,
        rc=rc,
        g=g,
        model="nonlinear",
        location="entry",
        modify_ref_part=modify_ref_part,
    )

    line = [dipedge1]

    # assign a segment
    sim.lattice.extend(line)

    # check the number of particles
    assert beam.total_number_of_particles() == 1

    # run simulation
    sim.track_particles()

    # collect final beam moments
    rbc = beam.beam_moments()

    sim.finalize()

    # return beam moments and reference particle
    return rbc, ref


def test_dipedge_modify_ref_part_false():
    """
    This tests the application of longitudinal particle boundary conditions.
    """
    edge_angle = math.pi / 8.0
    rc = 10.0
    g = 1.0e-3

    K0 = np.pi**2 / 6.0
    c2 = (1.0 / np.cos(edge_angle)) ** 3 * np.sin(edge_angle) / 2.0 * g**2 / rc**2 * K0
    c3 = (1.0 / np.cos(edge_angle)) ** 2 * g**2 / rc * K0
    c5 = np.tan(edge_angle) / (2.0 * rc)
    c13 = 0.0

    x_shift = -c3
    px_shift = c13 - c2 - c3 * c5

    rbc, ref = _track_particle_dipedge(False, edge_angle, g, rc)

    # access centroid data (= single particle orbit)
    vec_part = np.array(
        [
            rbc["mean_x"],
            rbc["mean_px"],
            rbc["mean_y"],
            rbc["mean_py"],
            rbc["mean_t"],
            rbc["mean_pt"],
        ]
    )
    vec_part_pred = [x_shift, px_shift, 0.0, 0.0, 0.0, 0.0]

    # access reference particle
    vec_ref = np.array([ref.x, ref.px, ref.y, ref.py])
    vec_ref_pred = [0, 0, 0, 0]

    rtol = 0.0
    atol = 1.0e-12
    np.testing.assert_allclose(vec_part, vec_part_pred, atol=atol)
    np.testing.assert_allclose(vec_ref, vec_ref_pred, atol=atol)


def test_dipedge_modify_ref_part_true():
    """
    This tests the application of longitudinal particle boundary conditions.
    """
    edge_angle = math.pi / 8.0
    rc = 10.0
    g = 1.0e-3

    K0 = np.pi**2 / 6.0
    c2 = (1.0 / np.cos(edge_angle)) ** 3 * np.sin(edge_angle) / 2.0 * g**2 / rc**2 * K0
    c3 = (1.0 / np.cos(edge_angle)) ** 2 * g**2 / rc * K0
    c5 = np.tan(edge_angle) / (2.0 * rc)
    c13 = 0.0
    x_shift = -c3
    px_shift = c13 - c2 - c3 * c5

    rbc, ref = _track_particle_dipedge(True, edge_angle, g, rc)

    # access centroid data (= single particle orbit)
    vec_part = np.array(
        [
            rbc["mean_x"],
            rbc["mean_px"],
            rbc["mean_y"],
            rbc["mean_py"],
            rbc["mean_t"],
            rbc["mean_pt"],
        ]
    )
    vec_part_pred = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # access reference particle
    vec_ref = np.array([ref.x, ref.px, ref.y, ref.py])
    vec_ref_pred = [x_shift, px_shift, 0, 0]

    rtol = 0.0
    atol = 1.0e-12
    np.testing.assert_allclose(vec_part, vec_part_pred, atol=atol)
    np.testing.assert_allclose(vec_ref, vec_ref_pred, atol=atol)
