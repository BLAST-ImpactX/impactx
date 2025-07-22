#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
from scipy.constants import c, e, m_e

from impactx import ImpactX, distribution, elements, twiss

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = True

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
distr = distribution.Triangle(
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

# design the accelerator lattice
ns = 1  # number of slices per ds in the element

# add beam diagnostics
monitor0 = elements.BeamMonitor("monitor0", backend="h5")
monitor1 = elements.BeamMonitor("monitor", backend="h5")

# element names consistent with HTU_base_lattice_matched.lte Elegant input

drift1 = elements.ExactDrift(name="drift1", ds=0.046, nslice=ns)
quad1 = elements.ExactQuad(
    name="quad1", ds=0.02903, k=207.0, unit=1, mapsteps=5, nslice=ns
)
quad2 = elements.ExactQuad(
    name="quad2", ds=0.02890, k=-207.0, unit=1, mapsteps=5, nslice=ns
)

# set the lattice
sim.lattice.append(monitor0)
sim.lattice.append(drift1)
sim.lattice.append(quad1)
sim.lattice.append(monitor1)
sim.lattice.append(quad2)
sim.lattice.append(monitor1)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
