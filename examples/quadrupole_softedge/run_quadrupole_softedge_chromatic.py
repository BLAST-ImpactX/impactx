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

kin_energy_MeV = 2.0e3  # reference energy

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
scale_min = 0.001
scale_max = 5.0
scale_arr = np.linspace(scale_min, scale_max, 50)
print(scale_arr)
if amr.ParallelDescriptor.IOProcessor():
    dx_arr = xi * scale_arr
    dpx_arr = pxi * scale_arr
    dy_arr = yi * scale_arr
    dpy_arr = pyi * scale_arr
    dt_arr = ti * scale_arr
    dpt_arr = pti * scale_arr
    pc.add_n_particles(
        dx_arr, dy_arr, dt_arr, dpx_arr, dpy_arr, dpt_arr, qm_eev, bunch_charge=0.0
    )

pc_in = sim.particle_container()

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
    gscale=2.0,
    cos_coefficients=cos_coeffs,
    sin_coefficients=sin_coeffs,
    mapsteps=200,
    nslice=ns,
)

quad_lin = elements.SoftQuadrupole(
    name="quad1",
    ds=1.5,
    gscale=2.0,
    cos_coefficients=cos_coeffs,
    sin_coefficients=sin_coeffs,
    mapsteps=200,
    nslice=ns,
)

# assign a lattice segment
sim.lattice.extend([quad_lin])

# run simulation
sim.track_particles()

# collective: gathers the particles of all ranks on the I/O rank
io_rank = amr.ParallelDescriptor.IOProcessor()
root_rank = amr.ParallelDescriptor.IOProcessorNumber()
pc_lin_f = sim.particle_container()
beam_lin_f = pc_lin_f.to_df(local=False, root_rank=root_rank)
if io_rank:
    beam_lin_f = beam_lin_f.set_index("idcpu").sort_index()

# set up second simulation
sim.lattice.clear()
pc.clear(keep_mass=True, keep_charge=True)

ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

#  add test particles
if amr.ParallelDescriptor.IOProcessor():
    dx_arr = xi * scale_arr
    dpx_arr = pxi * scale_arr
    dy_arr = yi * scale_arr
    dpy_arr = pyi * scale_arr
    dt_arr = ti * scale_arr
    dpt_arr = pti * scale_arr
    pc.add_n_particles(
        dx_arr, dy_arr, dt_arr, dpx_arr, dpy_arr, dpt_arr, qm_eev, bunch_charge=0.0
    )

sim.lattice.extend([quad_chr])

# run simulation
sim.track_particles()

pc_chr_f = sim.particle_container()
beam_chr_f = pc_chr_f.to_df(local=False, root_rank=root_rank)
if io_rank:
    beam_chr_f = beam_chr_f.set_index("idcpu").sort_index()

# clean shutdown
sim.finalize()

# analyze on the I/O rank, which holds all particles after the gather
if io_rank:
    # The particle ids of the second run continue after those of the first:
    # compare the runs in id order (the order of scale_arr), not by id value.
    phase_space = [
        "position_x",
        "momentum_x",
        "position_y",
        "momentum_y",
        "position_t",
        "momentum_t",
    ]
    diff = beam_chr_f[phase_space].to_numpy() - beam_lin_f[phase_space].to_numpy()

    # norm of the difference in the final phase space 6-vector
    norm_phase_space_diff = np.sqrt(np.sum(diff**2, axis=1))

    m, b = np.polyfit(np.log(scale_arr), np.log(norm_phase_space_diff), deg=1)
    print(f"  slope={m}")

    # Test for  quadratic scaling with initial phase space vector:
    rtol = 2.0e-2
    print(f"  rtol={rtol}")
    assert np.allclose(
        [m],
        [
            2.0,
        ],
        rtol=rtol,
    )

# plt.plot(scale_arr, norm_phase_space_diff)

# plt.show()
