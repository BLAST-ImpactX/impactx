#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.spin = True
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# basic beam parameters
kin_energy_MeV = 20.00  # reference energy (kinetic)
bunch_charge_C = 25.0e-12  # used with space charge
npart = 100000  # number of macro particles

# set reference particle
ref = sim.beam.ref
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
sigmaX = 0.003
sigmaY = sigmaX
sigmaPx = 0.2
sigmaPy = sigmaPx
sigmaT = 0.0

polarization_x = 0.4
polarization_y = 0.8
polarization_z = 0.2

distr = distribution.Gaussian(
    lambdaX=sigmaX,
    lambdaY=sigmaY,
    lambdaT=sigmaT,
    lambdaPx=sigmaPx,
    lambdaPy=sigmaPy,
    lambdaPt=0.0,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
spin = distribution.SpinvMF(
    polarization_x,
    polarization_y,
    polarization_z,
)

sim.add_particles(bunch_charge_C, distr, npart, spin)

# design the accelerator lattice
ns = 1  # number of slices per ds in the element

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

beta_gamma = ref.beta_gamma

ez_value = 10.0
bz_value = 3.0
ds_value = 50.0

sol_chr = elements.ChrAcc(
    name="sol_chr", ds=ds_value, ez=ez_value, bz=bz_value, nslice=ns
)

# set the lattice
sim.lattice.append(monitor)
sim.lattice.append(sol_chr)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
