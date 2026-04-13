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
sim.space_charge = False
sim.spin = True
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# basic beam parameters
kin_energy_MeV = 2.0e3  # reference energy (kinetic)
bunch_charge_C = 1.0e-9  # used with space charge
npart = 100000  # number of macro particles

# set reference particle
ref = sim.particle_container().push_ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    lambdaX=5.0e-6,
    lambdaY=8.0e-6,
    lambdaT=0.0599584916,
    lambdaPx=2.5543422003e-9,
    lambdaPy=1.5964638752e-9,
    lambdaPt=9.0e-4,
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

# parameters
L = 1.5
kquad = 60.0
ksol = 1.0

quad1 = elements.Quad(name="quad1", ds=L, k=kquad, nslice=ns)
iquad1 = elements.Quad(name="iquad1", ds=-L, k=kquad, nslice=ns)
quad2 = elements.ChrQuad(name="quad2", ds=L, k=-kquad, nslice=ns)
iquad2 = elements.ChrQuad(name="iquad2", ds=-L, k=-kquad, nslice=ns)
sol1 = elements.Sol(name="sol1", ds=L, ks=ksol, nslice=ns)
isol1 = elements.Sol(name="isol1", ds=-L, ks=ksol, nslice=ns)

# set the lattice
sim.lattice.append(monitor)
sim.lattice.append(quad1)
sim.lattice.append(iquad1)
sim.lattice.append(sol1)
sim.lattice.append(isol1)
sim.lattice.append(quad2)
sim.lattice.append(iquad2)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
