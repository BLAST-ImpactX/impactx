#!/usr/bin/env python3
#
# Copyright 2022-2023 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import matplotlib.pyplot as plt
import numpy as np
from conftest import basepath

from impactx import Config, ImpactX, amr, flatten_charge_to_2D


def test_charge_deposition_2D(save_png=True):
    """
    Deposit charge and access/plot it
    """
    sim = ImpactX()

    sim.n_cell = [16, 24, 32]
    sim.particle_shape = 2  # B-spline order
    sim.load_inputs_file(basepath + "/examples/fodo/input_fodo.in")
    sim.slice_step_diagnostics = False

    sim.init_grids()
    assert sim.n_cell == [16, 24, 32]

    sim.init_beam_distribution_from_inputs()
    sim.init_lattice_elements_from_inputs()

    # sim.track_particles()

    sim.deposit_charge()
    # rho = sim.rho(lev=0)
    # rs = rho.sum_unique(comp=0, local=False)

    rho_2d = flatten_charge_to_2D(sim)
    rho_2d_lvl0 = rho_2d[0][1]  # level 0 and unique boxes only
    rs = rho_2d_lvl0.sum_unique(comp=0, local=False)

    gm = sim.Geom(lev=0)
    dr = gm.data().CellSize()
    dV = np.prod(dr)

    beam_charge = dV * rs  # in C
    if Config.precision == "SINGLE":
        assert math.isclose(beam_charge, -1.0e-9, rel_tol=1e-6)
    else:
        assert math.isclose(beam_charge, -1.0e-9, rel_tol=1e-8)

    # plot data slices
    f = plt.figure()
    ax = f.gca()

    # comp = 0
    # mu = 1.0e6  # m->mu
    im = ax.imshow(
        rho_2d_lvl0[()] * dV,
        origin="lower",
        aspect="auto",
        # extent=[rbx.lo(1) * mu, rbx.hi(1) * mu, rbx.lo(0) * mu, rbx.hi(0) * mu],
    )
    cb = f.colorbar(im)
    cb.set_label(r"projected charge density  [C/m$^3$]")
    # ax.set_xlabel(r"$x$  [$\mu$m]")
    # ax.set_ylabel(r"$y$  [$\mu$m]")
    if save_png:
        plt.savefig("charge_deposition.png")
    else:
        plt.show()

    # finalize simulation
    sim.finalize()


# implement a direct script run mode, so we can run this directly too,
# with interactive matplotlib windows, w/o pytest
if __name__ == "__main__":
    amr.initialize([])

    test_charge_deposition_2D(save_png=False)

    if amr.initialized():
        amr.finalize()
