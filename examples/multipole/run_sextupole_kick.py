#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Ryan Sandberg, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

import amrex.space3d as amr
from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of  nm
kin_energy_MeV = 2.0e3  # reference energy

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = ref.charge_qe / (ref.mass_MeV * 1.0e6)  # electron charge/mass in e / eV

# set test particles
pc = sim.particle_container()

#  add test particles
xmin = 0.0
xmax = 5.0
if amr.ParallelDescriptor.IOProcessor():
    dx = np.linspace(xmin, xmax, 30)
    zero_arr = np.linspace(0, 0.0, 30)
    pc.add_n_particles(
        dx, zero_arr, zero_arr, zero_arr, zero_arr, zero_arr, qm_eev, bunch_charge=0.0
    )

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
multipole = [
    monitor,
    elements.Multipole(name="thin_sextupole", multipole=3, K_normal=10.0, K_skew=-2.0),
    monitor,
]
# assign a fodo segment
sim.lattice.extend(multipole)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
