#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

# from elements import LinearTransport

from impactx import ImpactX, Map3x6, Vector3, distribution, elements, twiss

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.spin = True
sim.slice_step_diagnostics = False

# domain decomposition & space charge mesh
sim.init_grids()

# set beam reference values
kin_energy_MeV = 10.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    **twiss(
        beta_x=2.8216194100262637,
        beta_y=2.8216194100262637,
        beta_t=0.5,
        emitt_x=2e-04,
        emitt_y=2e-04,
        emitt_t=2e-12,
        alpha_x=-1.5905003499999992,
        alpha_y=1.5905003499999992,
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

# design the accelerator lattice

# Elements for forward tracking
Q1 = elements.Quad(name="Q1", ds=1.0, k=1.0)
B1 = elements.Sbend(name="B1", ds=1.0, rc=10.0)

# Inverse spin map for Quad1
vmatQ1 = Vector3()
AmatQ1 = Map3x6()
AmatQ1[1, 3] = 27.846376647658047
AmatQ1[1, 4] = 12.868288416407257
AmatQ1[2, 1] = 19.9386437894808211
AmatQ1[2, 2] = 10.8925307463018033
Map1 = elements.SpinMap(v=vmatQ1, A=AmatQ1)

# Inverse spin map for Bend1
vmatB1 = Vector3()
AmatB1 = Map3x6()
vmatB1[2] = 2.269498669527234
AmatB1[1, 4] = 1.643140677773437
AmatB1[2, 1] = 0.236555147919017
AmatB1[2, 2] = 0.118376237268958
AmatB1[2, 6] = 0.096052815713841
AmatB1[3, 4] = -0.765638392807848
Map2 = elements.SpinMap(v=vmatB1, A=AmatB1)

line = [monitor, Map1, Q1, Map2, B1, monitor]

sim.lattice.extend(line)

# number of periods through the lattice
sim.periods = 1

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
