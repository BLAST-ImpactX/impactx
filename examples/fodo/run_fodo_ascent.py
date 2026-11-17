#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements, insitu

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = False

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

#   particle bunch
distr = distribution.Waterbag(
    lambdaX=3.9984884770e-5,
    lambdaY=3.9984884770e-5,
    lambdaT=1.0e-3,
    lambdaPx=2.6623538760e-5,
    lambdaPy=2.6623538760e-5,
    lambdaPt=2.0e-3,
    muxpx=-0.846574929020762,
    muypy=0.846574929020762,
    mutpt=0.0,
)
sim.add_particles(bunch_charge_C, distr, npart)

# in situ visualization: the phase space after each element
anim = insitu.Animation("phasespace", num_bins=64)
sim.hook["after_element"] = anim

# design the accelerator lattice
ns = 25  # number of slices per ds in the element
fodo = [
    elements.Drift(name="drift1", ds=0.25, nslice=ns),
    elements.Quad(name="quad1", ds=1.0, k=1.0, nslice=ns),
    elements.Drift(name="drift2", ds=0.5, nslice=ns),
    elements.Quad(name="quad2", ds=1.0, k=-1.0, nslice=ns),
    elements.Drift(name="drift3", ds=0.25, nslice=ns),
]
# assign a fodo segment
sim.lattice.extend(fodo)

# run simulation
sim.track_particles()

# the beam at the end, in 3D
fig = sim.beam.plot_ascent("beam3d", color="momentum_t", image_size=(800, 600))

# on MPI rank 0: save the images
if fig is not None:
    fig.savefig("diags/fodo_ascent_beam3d.png")
anim.save("diags/fodo_ascent_frames")  # one PNG per frame and projection
anim.save("diags/fodo_ascent.html")  # an interactive animation

# clean shutdown
sim.finalize()
