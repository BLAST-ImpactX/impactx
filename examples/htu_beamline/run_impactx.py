#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# HTU References:
# - https://doi.org/10.1103/vh62-gz1p
# - https://doi.org/10.1117/12.3056776
#
# -*- coding: utf-8 -*-

import numpy as np
from htu_lattice import get_lattice
from scipy.constants import c, e, m_e

from impactx import ImpactX, distribution, twiss

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = True

# silent running
silent = False
if silent:
    sim.verbose = 0
    sim.tiny_profiler = False
    # note: lattice beam monitors will still write files
    sim.diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# basic beam parameters
total_energy_MeV = 100.0  # reference energy (total)
mass_MeV = 0.510998950  # particle mass
kin_energy_MeV = total_energy_MeV - mass_MeV
bunch_charge_C = 25.0e-12  # used with space charge
npart = 10000  # number of macro particles

# set reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(mass_MeV).set_kin_energy_MeV(kin_energy_MeV)

# factors converting the beam distribution to ImpactX input
gamma = total_energy_MeV / mass_MeV
bg = np.sqrt(gamma**2 - 1.0)
sigma_tau = 1e-6  # in m
sigma_p = 2.5e-2  # dimensionless
rigidity = m_e * c * bg / e

#   particle bunch
distr = distribution.Gaussian(
    **twiss(
        beta_x=0.002,
        beta_y=0.002,
        beta_t=sigma_tau / sigma_p,
        emitt_x=1.5e-6 / bg,
        emitt_y=1.5e-6 / bg,
        emitt_t=sigma_tau * sigma_p,
        alpha_x=0.0,
        alpha_y=0.0,
        alpha_t=0.0,
    )
)
sim.add_particles(bunch_charge_C, distr, npart)

# set the lattice
sim.lattice.extend(get_lattice("impactx"))

# run simulation
sim.track_particles()

# in situ calculate the reduced beam characteristics
# note: rbc us a dataframe, sim can be finalized once rbc was created
# beam = sim.particle_container()
# rbc = beam.reduced_beam_characteristics()

# clean shutdown
sim.finalize()
