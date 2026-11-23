#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Nikita Kuklev
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
This test checks the physical invariant that survives the solve: the charge
density ``rho`` left on the grid after a space-charge solve must still integrate
to the beam charge. Using the integrated charge (``dV * sum(rho)``)
makes the check independent of the mesh resolution, which changes between solves.
"""

import math

import numpy as np
import pytest

from impactx import Config, ImpactX, distribution, elements

# vacuum permittivity [F/m], matching ablastr::constant::SI::epsilon_0.
EPSILON_0 = 8.8541878188e-12

BUNCH_CHARGE_C = 1.0e-8


def _build_sim(poisson_solver, prob_relative):
    """Set up a short constant-focusing channel with 3D space charge."""
    sim = ImpactX()

    sim.n_cell = [48, 48, 40]
    sim.tiny_profiler = False
    sim.particle_shape = 2
    sim.space_charge = "3D"
    sim.poisson_solver = poisson_solver
    sim.prob_relative = [prob_relative]
    sim.slice_step_diagnostics = False

    sim.init_grids()

    # 2 GeV proton beam (matches the cfchannel example / test_space_charge_fields)
    ref = sim.beam.ref
    ref.set_species("proton").set_kin_energy_MeV(2.0e3)

    distr = distribution.Waterbag(
        lambdaX=1.2154443728379865788e-3,
        lambdaY=1.2154443728379865788e-3,
        lambdaT=4.0956844276541331005e-4,
        lambdaPx=8.2274435782286157175e-4,
        lambdaPy=8.2274435782286157175e-4,
        lambdaPt=2.4415943602685364584e-3,
    )
    sim.add_particles(BUNCH_CHARGE_C, distr, 10000)

    # a single short CF cell, one slice -> exactly one space-charge solve
    sim.lattice.extend(
        [
            elements.ConstF(name="cf1", ds=2.0, kx=1.0, ky=1.0, kt=1.0, nslice=1),
        ]
    )
    return sim


def _integrated_charge(sim):
    """Return dV * sum(rho) on level 0, i.e. the total charge on the grid [C]."""
    rho = sim.rho(lev=0)
    rho_sum = rho.sum_unique(comp=0, local=False)  # MPI-collective
    cell_size = sim.Geom(lev=0).data().CellSize()
    dV = float(np.prod(cell_size))
    return dV * rho_sum


@pytest.mark.parametrize(
    "poisson_solver,prob_relative",
    [
        ("fft", 1.1),
        ("multigrid", 3.0),
    ],
)
def test_rho_integrates_to_beam_charge(poisson_solver, prob_relative):
    """rho after a SC solve must still integrate to the beam charge."""
    sim = _build_sim(poisson_solver, prob_relative)

    # one space-charge slice: this runs PoissonSolve on the deposited rho.
    # The "fft" (IGF) solver is only available when ImpactX was built with FFT
    # support; skip that parametrization on builds without it.
    try:
        sim.track_particles()
    except Exception as e:
        if poisson_solver == "fft" and "fft" in str(e).lower():
            sim.finalize()
            pytest.skip("ImpactX built without FFT support; skipping IGF solver")
        raise

    # integrated charge of the rho left on the grid by the space-charge solve
    q_on_grid = _integrated_charge(sim)

    # computePhi leaves rho alone, so the deposited beam charge is still on the
    # grid (up to the solver's deposition accuracy).
    rel_tol = 1e-4 if Config.precision == "SINGLE" else 1e-6
    assert math.isclose(abs(q_on_grid), BUNCH_CHARGE_C, rel_tol=rel_tol), (
        f"rho does not integrate to the beam charge after PoissonSolve "
        f"({poisson_solver}): integrated charge = {q_on_grid!r} C, expected "
        f"magnitude {BUNCH_CHARGE_C:.3e} C. A value near "
        f"{BUNCH_CHARGE_C * EPSILON_0:.3e} C would indicate rho was scaled "
        f"by epsilon_0."
    )

    sim.finalize()


if __name__ == "__main__":
    test_rho_integrates_to_beam_charge("fft", 1.1)
    test_rho_integrates_to_beam_charge("multigrid", 3.0)
