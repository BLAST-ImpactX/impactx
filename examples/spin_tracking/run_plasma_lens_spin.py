#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

from impactx import ImpactX, distribution, elements, twiss

# Physical reference parameters
APL_length = 20e-3  # [m]
APL_g = -1000
kin_energy_MeV = 200  # [MeV] reference energy
mass_MeV = 0.510998950  # [MeV]
bunch_charge_C = 1.0e-9  # used with space charge

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.spin = True
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

## Physics parameters for test (APL_g from input arguments)

# Load a 200 MeV electron beam with alpha=0 (x and y)
#  in the center of the APL
# reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#   particle bunch
distr = distribution.Gaussian(
    **twiss(
        beta_x=0.392,
        beta_y=0.392,
        beta_t=1.0,
        emitt_x=2.56e-8,
        emitt_y=2.56e-8,
        emitt_t=1e-06,
        alpha_x=0.02,
        alpha_y=0.02,
        alpha_t=0.0,
    )
)
spin_vectors = distribution.SpinvMF(
    0.6,
    0.5,
    0.4,
)
npart = 100000  # number of macro particles
sim.add_particles(bunch_charge_C, distr, npart, spin_vectors)

beam_in = sim.beam
rbc_in = beam_in.beam_moments()

# create the accelerator lattice

ns = 1  # number of slices per ds in the element
monitor = elements.BeamMonitor("monitor", backend="h5")

APL = elements.ChrPlasmaLens(name="APL", ds=APL_length, k=APL_g, unit=1, nslice=ns)

num_kicks = 20
APL_length_slice = APL_length / num_kicks

ThinTPL = elements.TaperedPL(name="TPL", k=-APL_g * APL_length_slice, taper=0.0, unit=1)
dr1 = elements.ChrDrift(name="dr", ds=-APL_length / (2 * num_kicks), nslice=ns)

invTPL = ([dr1, ThinTPL, dr1]) * (num_kicks)

# Set the lattice
sim.lattice.append(monitor)
sim.lattice.append(APL)
#    sim.lattice.extend(invTPL)
sim.lattice.append(monitor)

# run simulation
sim.track_particles()

# in situ calculate the reduced beam characteristics
beam_out = sim.beam
rbc_out = beam_out.beam_moments()

# clean shutdown
sim.finalize()

# access particle data
sigmaxi = rbc_in["sigma_x"]
sigmayi = rbc_in["sigma_y"]
sigmati = rbc_in["sigma_t"]
emittancexi = rbc_in["emittance_x"]
emittanceyi = rbc_in["emittance_y"]
emittanceti = rbc_in["emittance_t"]
meansxi = rbc_in["mean_sx"]
meansyi = rbc_in["mean_sy"]
meanszi = rbc_in["mean_sz"]

sigmaxf = rbc_out["sigma_x"]
sigmayf = rbc_out["sigma_y"]
sigmatf = rbc_out["sigma_t"]
emittancexf = rbc_out["emittance_x"]
emittanceyf = rbc_out["emittance_y"]
emittancetf = rbc_out["emittance_t"]
meansxf = rbc_out["mean_sx"]
meansyf = rbc_out["mean_sy"]
meanszf = rbc_out["mean_sz"]

# check that all final particle lie within the bucket
#    np.testing.assert_allclose(
#        [sigmaxf, sigmayf, sigmatf, emittancexf, emittanceyf, emittancetf],
#        [sigmaxi, sigmayi, sigmati, emittancexi, emittanceyi, emittanceti],
#        atol=1.0e-8,
#        rtol=0,
#    )
#    np.testing.assert_allclose(
#        [meansxi, meansyi, meanszi],
#        [ 0.6, 0.5, 0.4 ],
#        atol=1.0e-9,
#        rtol=0,
#    )
np.testing.assert_allclose(
    [meansxf, meansyf, meanszf],
    [0.6, 0.5, 0.4],
    atol=1.0e-9,
    rtol=0,
)
