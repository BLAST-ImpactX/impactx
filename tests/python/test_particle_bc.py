#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, distribution, elements


def test_particle_bc():
    """
    This tests the application of longitudinal particle boundary conditions.
    """
    sim = ImpactX()

    sim.particle_shape = 2
    sim.slice_step_diagnostics = False
    sim.particle_bc = "periodic"
    sim.init_grids()

    # init particle beam
    kin_energy_MeV = 10.0e3
    bunch_charge_C = 1.0e-9
    npart = 10000

    #   reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

    #   set bunch bucket length
    sim.particle_container().set_bucket_length(0.2)

    #   particle bunch
    distr = distribution.Gaussian(
        lambdaX=5.0e-4,
        lambdaY=5.0e-4,
        lambdaT=0.2,
        lambdaPx=0.0,
        lambdaPy=0.0,
        lambdaPt=0.0,
        muxpx=0.0,
        muypy=0.0,
        mutpt=0.0,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    pc = sim.particle_container()
    assert pc.total_number_of_particles() == npart

    # init accelerator lattice
    monitor = elements.BeamMonitor("monitor", backend="h5")

    lattice = [monitor, monitor]
    sim.lattice.extend(lattice)

    sim.track_particles()

    # access particle data
    [xmin, ymin, zmin, xmax, ymax, zmax] = pc.min_and_max_positions()

    # finalize simulation
    sim.finalize()

    # numerical parameters
    bucket_length = 0.2
    half_bucket_length = bucket_length / 2.0

    # check that all final particle lie within the bucket
    np.testing.assert_array_less([zmax], [half_bucket_length])
