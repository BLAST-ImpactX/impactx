#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np

import amrex.space3d as amr
from impactx import ImpactX, elements, fourier_coefficients

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True

# domain decomposition & space charge mesh
sim.init_grids()

# load a 2 GeV electron beam with an initial
# unnormalized rms emittance of 2 nm
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.beam.ref
ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)
qm_eev = ref.charge_qe / (ref.mass_MeV * 1.0e6)  # electron charge/mass in e / eV

#  set 6-vector of initial conditions
xi = 0.001
pxi = -0.01
yi = 0.02
pyi = 0.002
ti = 0.003
pti = -0.03

# set test particles
pc = sim.particle_container()

#  add test particles
scale_min = 0.0
scale_max = 10
if amr.ParallelDescriptor.IOProcessor():
    scale = np.linspace(scale_min, scale_max, 50)
    dx_arr = xi * scale
    dpx_arr = pxi * scale
    dy_arr = yi * scale
    dpy_arr = pyi * scale
    dt_arr = ti * scale
    dpt_arr = pti * scale
    pc.add_n_particles(
        dx_arr, dy_arr, dt_arr, dpx_arr, dpy_arr, dpt_arr, qm_eev, bunch_charge=0.0
    )

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 10  # number of slices per ds in the element

# specify the on-axis field profile
zmin = -0.75  # lower value of on-axis longitudinal coordinate (in meters)
zmax = 0.75  # upper value of on-axis longitudinal coordinate (in meters)
nz = 401  # number of longitudinal sampling points to be used
g = 0.07  # gap parameter (in meters)
L_hardedge = 1.0  # length in the hard-edge limit
zdata = np.linspace(zmin, zmax, nz)
bdata = (
    1.0
    / 2.0
    * (
        np.tanh((zdata + L_hardedge / 2.0) / g)
        - np.tanh((zdata - L_hardedge / 2.0) / g)
    )
)

# determine Fourier coefficients
ncoef = 25
cos_coeffs, sin_coeffs = fourier_coefficients(zdata, bdata, ncoef)

# design the accelerator lattice
quad_chr = elements.ChrSoftQuadrupole(
    name="quad1",
    ds=1.5,
    gscale=1.0,
    cos_coefficients=cos_coeffs,
    sin_coefficients=sin_coeffs,
    mapsteps=200,
    nslice=ns,
)

quad_lin = elements.SoftQuadrupole(
    name="quad1",
    ds=1.5,
    gscale=1.0,
    cos_coefficients=cos_coeffs,
    sin_coefficients=sin_coeffs,
    mapsteps=200,
    nslice=ns,
)

quad_chr_inv = elements.ChrSoftQuadrupole(
    name="quad1",
    ds=-1.5,
    gscale=1.0,
    cos_coefficients=cos_coeffs,
    sin_coefficients=sin_coeffs,
    mapsteps=200,
    nslice=ns,
)

quad_lin_inv = elements.SoftQuadrupole(
    name="quad1",
    ds=-1.5,
    gscale=1.0,
    cos_coefficients=cos_coeffs,
    sin_coefficients=sin_coeffs,
    mapsteps=200,
    nslice=ns,
)


# assign a fodo segment
sim.lattice.extend([monitor, quad_chr, monitor])

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
