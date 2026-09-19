#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, elements


def test_cfbend_zero_quad():
    """This test compares the linear map for an ExactCFbend element with zero quad strength to the linear map for an Sbend."""

    sim = ImpactX()

    # set numerical parameters and IO control
    sim.space_charge = False
    sim.slice_step_diagnostics = False

    # domain decomposition & space charge mesh
    sim.init_grids()

    #  set reference particle
    ref = sim.beam.ref
    kin_energy_MeV = 2.606299137493995  # reference kinetic energy (p = 70 MeV/c)
    ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

    # init accelerator lattice
    ns = 10  # number of slices per ds in the element

    # Construct lattice
    cfbend = elements.ExactCFbend(
        aperture_x=0.0,
        aperture_y=0.0,
        ds=0.391140372489,
        dx=0.0,
        dy=0.0,
        int_order=2,
        k_normal=[1.3386467172031062, 0.0, 30.21694, -75.81527],
        k_skew=[0.0, 0.0, 0.0, 0.0],
        mapsteps=10,
        name="m1r",
        nslice=ns,
        rotation=0.0,
        unit=0,
    )
    sbend = elements.Sbend(ds=0.391140372489, rc=0.747023084693581)

    sim.lattice.extend([cfbend])

    # return the linear map
    R_cfbend = sim.lattice.transfer_map(ref).to_numpy()
    np.savetxt("matrix_cfbend.txt", R_cfbend, delimiter=" ")

    sim.lattice.clear()
    sim.lattice.extend([sbend])

    # return the linear map
    R_sbend = sim.lattice.transfer_map(ref).to_numpy()
    np.savetxt("matrix_sbend.txt", R_cfbend, delimiter=" ")

    atol = 1.0e-14
    rtol = 0.0
    print(f"  atol={atol} (ignored: rtol~={rtol})")
    print("  Linear map for cfbend with zero quad strength: ")
    print(R_cfbend)
    print()
    print("  Linear map for sbend:  ")
    print(R_sbend)

    assert np.allclose(
        R_cfbend,
        R_sbend,
        rtol=rtol,
        atol=atol,
    )

    # clean shutdown
    sim.finalize()
