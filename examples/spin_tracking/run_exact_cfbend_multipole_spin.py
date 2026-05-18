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
sim.spin = True
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# basic beam parameters
kin_energy_MeV = 2.0e3
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

# set reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

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
spin = distribution.SpinvMF(
    0.4,
    0.9,
    0.1,
)

sim.add_particles(bunch_charge_C, distr, npart, spin)

# design the accelerator lattice
ns = 1  # number of slices per ds in the element

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# lattice elements
cfbend1 = elements.ExactCFbend(
    name="cfbend1",
    ds=1.0,
    k_normal=[1.0e-7, 1.0, -2.0],
    k_skew=[0.0, -0.5, 1.4],
    unit=0,
    int_order=2,
    mapsteps=20,
    nslice=ns,
)

cfbend2 = elements.ExactCFbend(
    name="cfbend1",
    ds=-1.0,
    k_normal=[1.0e-7, 1.0, -2.0],
    k_skew=[0.0, -0.5, 1.4],
    unit=0,
    int_order=2,
    mapsteps=20,
    nslice=ns,
)

multipole1 = elements.ExactMultipole(
    name="quad1",
    ds=-1.0,
    k_normal=[0.0, 1.0, -2.0],
    k_skew=[0.0, -0.5, 1.4],
    unit=0,
    int_order=2,
    mapsteps=20,
    nslice=ns,
)

multipole2 = elements.ExactMultipole(
    name="quad1",
    ds=1.0,
    k_normal=[0.0, 1.0, -2.0],
    k_skew=[0.0, -0.5, 1.4],
    unit=0,
    int_order=2,
    mapsteps=20,
    nslice=ns,
)

# set the lattice
sim.lattice.append(monitor)
sim.lattice.append(cfbend1)
# sim.lattice.append(cfbend2)
sim.lattice.append(multipole1)
# sim.lattice.append(multipole2)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
