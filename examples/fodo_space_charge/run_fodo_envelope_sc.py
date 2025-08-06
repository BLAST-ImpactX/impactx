#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements, twiss

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 0  # B-spline order
sim.space_charge = "2D"
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# model a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 6.7  # reference energy

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(938.27208816).set_kin_energy_MeV(kin_energy_MeV)

#  beam current in A
beam_current_A = 0.5

#   particle bunch
distr = distribution.KVdist(
    **twiss(
        beta_x=0.737881,
        beta_y=0.737881,
        beta_t=0.5,
        emitt_x=1.0e-6,
        emitt_y=1.0e-6,
        emitt_t=1.0e-12,
        alpha_x=2.4685083,
        alpha_y=-2.4685083,
        alpha_t=0.0,
    )
)

sim.init_envelope(ref, distr, beam_current_A)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 50  # number of slices per ds in the element
fodo = [
    monitor,
    elements.Drift(name="drift1", ds=7.44e-2, nslice=ns),
    elements.Quad(name="quad1", ds=6.10e-2, k=-103.12574100336, nslice=ns),
    elements.Drift(name="drift2", ds=14.88e-2, nslice=ns),
    elements.Quad(name="quad2", ds=6.10e-2, k=103.12574100336, nslice=ns),
    elements.Drift(name="drift3", ds=7.44e-2, nslice=ns),
    monitor,
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.track_envelope()

# clean shutdown
sim.finalize()
