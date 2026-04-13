#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 230 MeV electron beam with an initial
# unnormalized rms emittance of 1 mm-mrad in all
# three phase planes
kin_energy_MeV = 230.0  # reference energy

#   reference particle
ref = sim.particle_container().push_ref_particle()
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

# design the accelerator lattice

# access RF cavity on-axis field data
data_in = np.loadtxt("onaxis_data.in")
z = data_in[:, 0]
ez_onaxis = data_in[:, 1]
ncoef = 25

#   Drift elements
dr1 = elements.Drift(name="dr1", ds=0.4, nslice=1)
dr2 = elements.Drift(name="dr2", ds=0.032997, nslice=1)
#   RF cavity element
rf = elements.RFCavity(
    name="rf",
    ds=1.31879807,
    escale=62.0,
    z=z,
    field_on_axis=ez_onaxis,
    ncoef=ncoef,
    freq=1.3e9,
    phase=85.5,
    mapsteps=100,
    nslice=4,
)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

sim.lattice.extend(
    [
        monitor,
        dr1,
        dr2,
        rf,
        dr2,
        dr2,
        rf,
        dr2,
        dr2,
        rf,
        dr2,
        dr2,
        rf,
        dr2,
        monitor,
    ]
)

# run simulation
sim.track_reference(ref)

# clean shutdown
sim.finalize()
