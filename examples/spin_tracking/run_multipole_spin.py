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
# sim.diagnostics = False  # benchmarking
sim.spin = True
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=3.9984884770e-4,
    lambdaY=3.9984884770e-4,
    lambdaT=1.0e-3,
    lambdaPx=2.6623538760e-4,
    lambdaPy=2.6623538760e-4,
    lambdaPt=2.0e-3,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
spin = distribution.SpinvMF(
    0.4,
    0.9,
    0.1,
)
sim.add_particles(bunch_charge_C, distr, npart, spin)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice)
ns = 5  # number of slices per ds in the element
Kn = 2.0
L = 1.0
Kn_integrated = 1.0
Ks_integrated = -0.5
L2 = 1.0e-4

lattice = [
    monitor,
    elements.ExactMultipole(
        name="multipole",
        ds=L,
        k_normal=[0.0, Kn],
        k_skew=[0.0, 0.0],
        mapsteps=100,
        nslice=ns,
    ),
    elements.ExactQuad(
        name="quad",
        ds=-L,
        k=Kn,
        mapsteps=100,
        nslice=ns,
    ),
    #  Test 2:  thick multipole vs thin multipole
    elements.Multipole(
        name="thin_multipole",
        multipole=3,
        K_normal=-Kn_integrated,
        K_skew=-Ks_integrated,
    ),
    elements.ExactMultipole(
        name="thick_multipole",
        ds=L2,
        k_normal=[0.0, 0.0, Kn_integrated / L2],
        k_skew=[0.0, 0.0, Ks_integrated / L2],
        mapsteps=100,
        nslice=ns,
    ),
    monitor,
]
sim.lattice.extend(lattice)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
