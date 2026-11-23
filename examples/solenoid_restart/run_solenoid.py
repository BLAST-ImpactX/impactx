#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import Config, ImpactX, elements, push

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load the particle bunch and the reference particle (a 250 MeV proton beam)
# from the final beam monitor output of the solenoid example
beam = sim.beam
push(beam, elements.Source("openPMD", "../solenoid.py/diags/openPMD/m1.h5"))

# test only:
# check that the reference particle was restored exactly from the file metadata
# this is a 250 MeV proton at the end of the first solenoid channel (one period)
rtol = 1.0e-13 if Config.precision == "DOUBLE" else 1.0e-5  # a few machine epsilon
ref = beam.ref
assert abs(ref.kin_energy_MeV - 250.0) < 250.0 * rtol
assert abs(ref.charge_qe - 1.0) < rtol
assert abs(ref.s - 3.820395) < 3.820395 * rtol

# add beam diagnostics
m1 = elements.BeamMonitor("m1", backend="h5")

# design the accelerator lattice
sol1 = elements.Sol(name="sol1", ds=3.820395, ks=0.8223219329893234)
sim.lattice.append(m1)
sim.lattice.extend([sol1] * 3)
sim.lattice.append(m1)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
