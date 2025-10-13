#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Chad Mitchell, Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from impactx import ImpactX, distribution, elements

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
sim.slice_step_diagnostics = True
sim.tiny_profiler = False

# domain decomposition & space charge mesh
sim.init_grids()

# load a 5 GeV electron beam with an initial
# normalized transverse rms emittance of 1 um
kin_energy_MeV = 2.0e3  # reference energy
bunch_charge_C = 1.0e-9  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

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
sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")

# design the accelerator lattice
ns = 1  # number of slices per ds in the element

bend = [
    monitor,
    elements.Drift(ds=1.0, nslice=ns, name="drift1"),
]

# assign a lattice segment
sim.lattice.extend(bend)
sim.periods = 10


# at the beginning of each period (e.g, turn or channel period), modify the lattice
def hook_before_period(sim):
    turn = sim.tracking_period
    step = sim.tracking_step
    print("- Before turn:")
    print(f"  Updating lattice at turn {turn}, step {step}", flush=True)

    beam = sim.particle_container()
    ref = beam.ref_particle()
    rbc = beam.beam_moments()
    print(
        f"  Beam at s={ref.s:.2f}m, t={ref.t:.2f}s with beta_x={rbc['beta_x']}m",
        flush=True,
    )

    if turn > 0:
        next_ds = 1.0 / ref.s
        print(f"  Updating next drift ds to: {next_ds:.2f}m", flush=True)
        sim.lattice.select(name="drift1")[0].ds = next_ds


def hook_before_element(sim):
    element = sim.tracking_element
    print("- Before element:")
    print(
        f"  Current element name: {element.name} with ds={element.ds:.2f}", flush=True
    )

    if element.name != "monitor":
        print(f"  Drift ds is: {element.ds:.2f}m", flush=True)


sim.hook["before_period"] = hook_before_period
sim.hook["before_element"] = hook_before_element

# run simulation
sim.track_particles()

# clean shutdown
sim.finalize()
