#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
from scipy import constants as sconst

from impactx import ImpactX, elements, fourier_coefficients

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# reference kinetic energy (initial)
kin_energy_MeV = 230.0  # reference energy

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

# design the accelerator lattice

# access RF cavity on-axis field data
data_in = np.loadtxt("onaxis_data.in")
z = data_in[:, 0]
ez_onaxis = data_in[:, 1]
ncoef = 25

cos_coeffs, sin_coeffs = fourier_coefficients(z, ez_onaxis, ncoef)

def change_in_gamma(cos_coeffs, sin_coeffs, f, L, emax, beta, phase, t0):
    # predicted energy gain
    zmin = -L/(2.0*beta) 
    theta = phase * np.pi / 180.0
    w = 2.0*np.pi*f
    k = w / sconst.c
    sinkL = np.sin(k*L/(2.0*beta))

    Abasearr = np.array([(-1)**j/((k*L-2*j*np.pi*beta)*(k*L+2*j*np.pi*beta)) for j in range(ncoef)])
    Bbasearr = np.array([(-1)**j*j/((k*L-2*j*np.pi*beta)*(k*L+2*j*np.pi*beta)) for j in range(ncoef)])
    Acoeffs = -2.0 * Abasearr * k * L**2 * beta * cos_coeffs * sinkL
    Bcoeffs = 4.0 * Bbasearr * np.pi * L * beta**2 * sin_coeffs * sinkL
    Asum = sum(Acoeffs[1:ncoef-1]) + Acoeffs[0]/2.0
    Bsum = sum(Bcoeffs[1:ncoef-1])

    dpt = emax * (Asum * np.cos(k*(t0-zmin/beta) + theta) + Bsum * np.sin(k*(t0-zmin/beta) + theta))
    dgamma = -dpt
    return dgamma

#   Drift elements
dr1 = elements.Drift(name="dr1", ds=0.4, nslice=1)
dr2 = elements.Drift(name="dr2", ds=0.032997, nslice=1)

#   RF cavity element
rf = elements.RFCavity(
    name="rf",
    ds=1.31879807,
    escale=62.0,
    z=z,
    field_or_gradient=ez_onaxis,
    ncoef=ncoef,
    freq=1.3e9, 
    phase=85.5,
    mapsteps=100,
    nslice=4,
)


# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

#sim.lattice.extend(
#    [
#        monitor,
#        dr1,
#        dr2,
#        rf,
#        dr2,
#        dr2,
#        monitor,
#    ]
#)

#sim.lattice.extend(
#    [
#        monitor,
#        rf,
#        monitor,
#    ]
#)

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

def hook_before_element(sim):
    element = sim.tracking_element
    if element.name == "rf":
        beam = sim.particle_container()
        ref = beam.ref_particle()
        beta = ref.beta
        f = element.freq
        L = element.ds
        emax = element.escale
        phase = element.phase
        t0 = ref.t
        dgamma = change_in_gamma(cos_coeffs, sin_coeffs, f, L, emax, beta, phase, t0)
        print(
            f"  Beam at s={ref.s:.2f}m, t={ref.t:.2f}s, gammai={ref.gamma:.2f}, gamma change={dgamma:.2f}",
            flush=True,
        )    
 
sim.hook["before_element"] = hook_before_element
    
# run simulation
sim.track_reference(ref)

# clean shutdown
sim.finalize()
