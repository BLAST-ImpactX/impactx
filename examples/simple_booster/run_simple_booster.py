#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Eric G. Stern, Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from booster_impactx_lattice import get_lattice
from scipy.constants import c, eV, m_p, pi

from impactx import ImpactX, distribution, elements, twiss

mp_mev = 1.0e-6 * m_p * c**2 / eV
total_Booster_charge = 6.7e12  # PIP-II full Booster
active_buckets = 81  # 81 out of 84 buckets full
harmonic_number = 84

turns = 2

emit_x = 16.0e-6  # normalized 95% emit
emit_y = 16.0e-6
emit_eV_s = 0.1  # longitudinal emittance 97% eV-s

#
# these lattice functions are calculated with Synergia3
# from the sbbooster-cooked.madx file.
alpha_x = -1.298673960026007664e-02
beta_x = 3.373645362843065243e01
alpha_y = 6.089861210659328755e-03
beta_y = 5.252517912567207681e00

disp_x = 3.187407765856291153e00
disp_px = 1.136005067625678322e-03
# s betax alphax psix dispx dprimex betay alphay psiy dispy dprimey
# 4.742027519999999186e+02 3.373645362843065243e+01 -1.298673960026007664e-02 4.216017341852963085e+01 3.187407765856291153e+00 1.136005067625678322e-03 5.252517912567207681e+00 6.089861210659328755e-03 4.278849194102631515e+01 0.000000000000000000e+00 0.000000000000000000e+00

sim = ImpactX()

# set numerical parameters and IO control
sim.space_charge = False
# sim.diagnostics = False  # benchmarking
sim.slice_step_diagnostics = True


# domain decomposition & space charge mesh
sim.init_grids()

# set the bucket length so that particles can be phase wrapped
booster_lattice = get_lattice()
total_s = sum(element.ds for element in booster_lattice)
print(f"Read Booster lattice, {len(booster_lattice)} elements, len = {total_s}")

bucket_length = total_s/harmonic_number
print("bucket_length: ", bucket_length)
sim.particle_container().set_bucket_length(bucket_length)

sim.particle_bc = "periodic"

# load a 800 MeV proton beam

kin_energy_MeV = 800.0  # reference energy 800 MeV
bunch_charge_C = eV * total_Booster_charge / active_buckets  # used with space charge
npart = 10000  # number of macro particles

#   reference particle
ref = sim.particle_container().ref_particle()
ref.set_species("proton").set_kin_energy_MeV(kin_energy_MeV)

distr = distribution.Gaussian(
    **twiss(
        beta_x=beta_x,
        alpha_x=alpha_x,
        emitt_x=emit_x / (ref.beta_gamma * 6),  # normalized 95% emit -> geometric
        beta_y=beta_y,
        alpha_y=alpha_y,
        emitt_y=emit_y / (ref.beta_gamma * 6),  # normalized 95% emit -> geometric
        emitt_t=emit_eV_s
        * 1.0e-6
        * c
        / (ref.mass_MeV * ref.beta_gamma * 6 * pi),  # 97% emit eV-s -> RMS emit
        beta_t=1258.0,  # seems to go from 1158 close to the center to
        # 1258 at about 1.25m
        # dispersion from Synergia so it needs conversion from dp/p to dE/p
        dispersion_x=disp_x / ref.beta,
        dispersion_px=disp_px / ref.beta,
    )
)

sim.add_particles(bunch_charge_C, distr, npart)

# add beam diagnostics
monitor = elements.BeamMonitor("monitor", backend="h5")
zerodrift = elements.Drift(ds=0.0)
# initial 0 length drift is to force initial phase wrapping
# run a lattice with a single monitor to capture the initial beam.
sim.lattice.clear()
sim.lattice.append(zerodrift)
sim.lattice.append(monitor)
sim.periods = 1
sim.track_particles()

# clear the lattice and load the Booster lattice
sim.lattice.clear()
sim.lattice.extend(get_lattice())
sim.lattice.append(monitor)

# run simulation
sim.periods = turns
sim.track_particles()

# clean shutdown
sim.finalize()
