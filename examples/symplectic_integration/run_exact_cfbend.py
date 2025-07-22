#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.particle_shape = 2  # B-spline order
sim.space_charge = False
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# basic beam parameters
mass_MeV = 0.510998950  # particle mass
kin_energy_MeV = 2.0e3
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

# set reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(mass_MeV).set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=4.0e-5,
    lambdaY=4.0e-5,
    lambdaT=1.0e-3,
    lambdaPx=3.0e-5,
    lambdaPy=3.0e-5,
    lambdaPt=2.0e-4,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# design the accelerator lattice
ns = 1  # number of slices per ds in the element

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# lattice elements
cfbend1 = elements.ExactCFbend(
    name="cfbend1",
    ds=1.0,
    k_normal=[0.1, 0.0],
    k_skew=[0.0, 0.0],
    unit=0,
    int_order=2,
    mapsteps=5,
    nslice=ns,
)
sbend1 = elements.ExactSbend(name="sbend1", ds=-1.0, phi=-5.729577951308232, nslice=ns)

# set the lattice
sim.lattice.append(monitor)
sim.lattice.append(cfbend1)
sim.lattice.append(sbend1)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
