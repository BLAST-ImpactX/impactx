#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

# from elements import LinearTransport

from impactx import ImpactX, Vector3, distribution, elements, twiss

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.spin = True
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 45.6e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
gyromagnetic_anomaly = 0.00115965218062  # value for an electron
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(
    kin_energy_MeV
).set_gyromagnetic_anomaly(gyromagnetic_anomaly)

#   target beta functions (m)
beta_star_x = 0.15
beta_star_y = 0.8e-3
beta_star_t = 9.210526315789473

#   particle bunch
distr = distribution.Waterbag(
    **twiss(
        beta_x=beta_star_x,
        beta_y=beta_star_y,
        beta_t=beta_star_t,
        emitt_x=0.27e-09,
        emitt_y=1.0e-12,
        emitt_t=1.33e-06,
        alpha_x=0.0,
        alpha_y=0.0,
        alpha_t=0.0,
    )
)
spin = distribution.SpinvMF(
    0.4,
    0.9,
    0.1,
)
sim.add_particles(bunch_charge_C, distr, npart, spin)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# initialize the spin map generators
vmat = Vector3

# matrix elements for the horizontal plane
vmat = [1.0, 0.0, 0.0]

Amat[1, 1] = 0.642252653176584
Amat[1, 2] = 0.114973951021402
Amat[2, 1] = -5.109953378728999
Amat[2, 2] = 0.642252653176584
Amat[3, 1] = 0.5
Amat[3, 6] = -0.2

# design the accelerator lattice
map = [
    monitor,
    elements.SpinMap(v=vmat, A=Amat),
]

sim.lattice.extend(map)

# number of periods through the lattice
sim.periods = 1

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
