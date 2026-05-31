#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import numpy as np

from impactx import ImpactX, elements


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
    beam.add_n_particles(dx, dy, dt, dpx, dpy, dpt, qm_eev, bunch_charge_C)

    # design the accelerator lattice)
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

    ref_state = ref.copy()

    sim.finalize()

    # return beam moments and reference particle
    return rbc, ref_state


def test_dipedge_modify_ref_part_false():
    """
    This tests the application of longitudinal particle boundary conditions.
    """
    edge_angle = math.pi / 8.0
    rc = 10.0
    g = 1.0e-2

    K0 = np.pi**2 / 6.0
    c2 = (1.0 / np.cos(edge_angle)) ** 3 * np.sin(edge_angle) / 2.0 * g**2 / rc**2 * K0
    c3 = (1.0 / np.cos(edge_angle)) ** 2 * g**2 / rc * K0
    c5 = np.tan(edge_angle) / (2.0 * rc)
    c13 = 0.0

    rbc, ref = _track_particle_dipedge(False, edge_angle, g, rc)

    x_shift = -c3
    px_shift = c13 - c2 - c3 * c5
    t_shift = -c3 * px_shift / ref.beta

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
    vec_part_pred = [x_shift, px_shift, 0.0, 0.0, t_shift, 0.0]

    # access reference particle
    vec_ref = np.array([ref.x, ref.px, ref.y, ref.py, ref.t])
    vec_ref_pred = [0, 0, 0, 0, 0]
    vec_ref_on_shell = (
        (ref.px) ** 2 + (ref.py) ** 2 + (ref.pz) ** 2 - (ref.pt) ** 2 + 1.0
    )

    atol = 1.0e-12 if Config.precision != "SINGLE" else 1.0e-6
    np.testing.assert_allclose(vec_part, vec_part_pred, atol=atol)
    np.testing.assert_allclose(vec_ref, vec_ref_pred, atol=atol)
    np.testing.assert_allclose(vec_ref_on_shell, 0.0, atol=atol)


def test_dipedge_modify_ref_part_true():
    """
    This tests the application of longitudinal particle boundary conditions.
    """
    edge_angle = math.pi / 8.0
    rc = 10.0
    g = 1.0e-2

    K0 = np.pi**2 / 6.0
    c2 = (1.0 / np.cos(edge_angle)) ** 3 * np.sin(edge_angle) / 2.0 * g**2 / rc**2 * K0
    c3 = (1.0 / np.cos(edge_angle)) ** 2 * g**2 / rc * K0
    c5 = np.tan(edge_angle) / (2.0 * rc)
    c13 = 0.0

    rbc, ref = _track_particle_dipedge(True, edge_angle, g, rc)

    x_shift = -c3
    px_shift = c13 - c2 - c3 * c5
    t_shift = -c3 * px_shift / ref.beta

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
    vec_ref = np.array([ref.x, ref.px, ref.y, ref.py, ref.t])
    vec_ref_pred = [x_shift, px_shift, 0, 0, t_shift]
    vec_ref_on_shell = (
        (ref.px) ** 2 + (ref.py) ** 2 + (ref.pz) ** 2 - (ref.pt) ** 2 + 1.0
    )
    atol = 1.0e-12 if Config.precision != "SINGLE" else 1.0e-6
    np.testing.assert_allclose(vec_part, vec_part_pred, atol=atol)
    np.testing.assert_allclose(vec_ref, vec_ref_pred, atol=atol)
    np.testing.assert_allclose(vec_ref_on_shell, 0.0, atol=atol)
    np.testing.assert_allclose(vec_ref_on_shell, 0.0, atol=atol)
