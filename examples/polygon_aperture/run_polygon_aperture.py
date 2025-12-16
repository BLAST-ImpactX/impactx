#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from scipy.constants import c, eV, m_p

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 800.0  # reference energy 800 MeV proton
bunch_charge_C = 1.0e-9  # used with space charge
npart = 50000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(1.0).set_mass_MeV(1.0e-6 * m_p * c**2 / eV).set_kin_energy_MeV(
    kin_energy_MeV
)

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=2.0e-3,
    lambdaY=2.0e-3,
    lambdaT=0.4,
    lambdaPx=4.0e-4,
    lambdaPy=4.0e-4,
    lambdaPt=2.0e-3,
    muxpx=0.0,
    muypy=0.0,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

vertices_x = [
    float(u)
    for u in "0.5e-3 0.5e-3 -0.5e-3 -0.5e-3 -1.5e-3 -1.5e-3 -0.5e-3 -0.5e-3 0.5e-3 0.5e-3 1.5e-3 1.5e-3 0.5e-3".split()
]
vertices_y = [
    float(u)
    for u in "0.5e-3 1.5e-3 1.5e-3 0.5e-3 0.5e-3 -0.5e-3 -0.5e-3 -1.5e-3 -1.5e-3 -0.5e-3 -0.5e-3 0.5e-3 0.5e-3".split()
]
mr2 = 2 * 0.5e-3**2

aperture = elements.PolygonAperture(vertices_x, vertices_y, action="transmit")

print(aperture.to_dict())
assert aperture.min_radius2 == 0.0
aperture.min_radius2 = mr2
assert abs(aperture.min_radius2 / mr2 - 1.0) < 1.0e-15

# design the accelerator lattice)
ns = 1  # number of slices per ds in the element
channel = [
    monitor,
    aperture,
    monitor,
]

sim.lattice.extend(channel)

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
