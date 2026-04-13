#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
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

# load a 5 GeV electron beam with an initial
# normalized transverse rms emittance of 1 um
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().push_ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=5.0e-6,  # 5 um
    lambdaY=8.0e-6,  # 8 um
    lambdaT=0.0599584916,  # 200 ps
    lambdaPx=2.5543422003e-9,  # exn = 50 pm-rad
    lambdaPy=1.5964638752e-9,  # eyn = 50 pm-rad
    lambdaPt=9.0e-4,  # approximately dE/E
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

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 1  # number of slices per ds in the element

# element length (m)
L = 1.5

# large curvature radius value for the quad test
rclarge = 1.0e8

# small gradient value for the sbend test
ksmall = 1.0e-8

lattice = [
    monitor,
    elements.CFbend(name="cfbend1", ds=L, rc=rclarge, k=60.0, nslice=ns),
    elements.Quad(name="iquad1", ds=-L, k=60.0, nslice=ns),
    elements.CFbend(name="cfbend2", ds=L, rc=10.0, k=ksmall, nslice=ns),
    elements.Sbend(name="isbend1", ds=-L, rc=10.0, nslice=ns),
    monitor,
]

# assign a lattice segment
sim.lattice.extend(lattice)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
