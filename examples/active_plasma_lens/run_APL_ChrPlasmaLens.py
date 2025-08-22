#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


def run_APL_ChrPlasmaLens(APL_g: float, sigpt_0: float):
    "Run the ChrPlasmaLens simulation with the given APL gradient APL_g [T/m] and sigma_pt [-]"
    import math

    from impactx import ImpactX, distribution, elements

    sim = ImpactX()

    # set numerical parameters and IO control
    sim.space_charge = False
    # sim.diagnostics = False  # benchmarking
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    ## Physics parameters for test (APL_g from input arguments)
    APL_length = 20e-3  # [m]

    # Load a 200 MeV electron beam with alpha=0 (x and y)
    #  in the center of the APL
    kin_energy_MeV = 200  # reference energy
    bunch_charge_C = 1.0e-9  # used with space charge
    # reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(0.510998950).set_kin_energy_MeV(kin_energy_MeV)

    # Midpoint parameters
    alpha_mid = 0.0
    sigma_mid = 0.5e-3 / 2  # [m]
    emitn = 10e-6  # [m]
    emitg = emitn / ref.beta_gamma
    beta_mid = sigma_mid**2 / emitg
    gamma_mid = 1 / beta_mid  # [1/m]
    print(
        f"sigma_mid = {sigma_mid} [m], beta_mid = {beta_mid} [m], gamma_mid = {gamma_mid} [m]"
    )
    print(
        f"emitn = {emitn} [m], emitg = {emitg} [m], ref.beta_gamma = {ref.beta_gamma}"
    )
    print()

    # Back-propagate 1/2 lens length as in vacuum
    # (from symmetry point in the middle of the lens)
    beta_0 = beta_mid + (APL_length / 2) ** 2 / beta_mid
    alpha_0 = +APL_length / 2 / beta_mid
    gamma_0 = gamma_mid
    sigma_0 = math.sqrt(emitg * beta_0)
    sigmap_0 = math.sqrt(emitg * gamma_0)
    mu_0 = alpha_0 / math.sqrt(beta_0 * gamma_0)
    print(
        f"sigma_0 = {sigma_0} [m], beta_0 = {beta_0} [m], alpha_0 = {alpha_0}, sigmap_0 = {sigmap_0}"
    )
    print()

    # Longitudinal parameters (sigpt_0 [-] from input arguments)
    sigt_0 = 1e-3  # [m]
    emit_t = math.sqrt(sigt_0**2 * sigpt_0**2 - 0**2)
    print(f"sigt_0 = {sigt_0} [m], sigpt_0 = {sigpt_0} [-], emit_t = {emit_t}")

    print()

    #   particle bunch
    distr = distribution.Waterbag(
        lambdaX=sigma_0,
        lambdaY=sigma_0,
        lambdaT=sigt_0,
        lambdaPx=sigmap_0,
        lambdaPy=sigmap_0,
        lambdaPt=sigpt_0,
        muxpx=mu_0,
        muypy=mu_0,
        mutpt=0.0,
    )
    npart = 10000  # number of macro particles
    sim.add_particles(bunch_charge_C, distr, npart)

    # design the accelerator lattice

    ns = 25  # number of slices per ds in the element
    monitor = elements.BeamMonitor("monitor", backend="h5")
    APL = elements.ChrPlasmaLens(name="APL", ds=APL_length, k=APL_g, unit=1, nslice=ns)

    lattice = [
        monitor,
        APL,
        monitor,
    ]
    # assign a fodo segment
    sim.lattice.extend(lattice)

    # run simulation
    sim.track_particles()

    # clean shutdown
    sim.finalize()
