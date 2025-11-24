#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
import math

import numpy as np

# Physical reference parameters
APL_length = 20e-3  # [m]
kin_energy_MeV = 200  # [MeV] reference energy
mass_MeV = 0.510998950  # [MeV]
bunch_charge_C = 1.0e-9  # used with space charge


def run_APL_tracking(
    APL_g: float, sigpt_0: float, sigma_mid: float, lensType: str = "ChrPlasmaLens"
):
    """
    Run a plasma lens tracking simulation with the given APL gradient APL_g [T/m], sigma_pt [-], and sigma_mid [m].
    Can use lensType='ChrPlasmaLens' | 'ConstF' | 'ChrDrift' (expects APL_g = 0) | 'ChrQuad' (only horizontal plane valid)
    """

    from impactx import ImpactX, distribution, elements

    print("", flush=True)
    print(
        f"*** run_APL_tracking({APL_g}, {sigpt_0}, {sigma_mid}, '{lensType}') :",
        flush=True,
    )
    print("", flush=True)

    sim = ImpactX()

    # set numerical parameters and IO control
    sim.space_charge = False
    # sim.diagnostics = False  # benchmarking
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    ## Physics parameters for test (APL_g from input arguments)

    # Load a 200 MeV electron beam with alpha=0 (x and y)
    #  in the center of the APL
    # reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(mass_MeV).set_kin_energy_MeV(kin_energy_MeV)

    (emitg, beta_0, gamma_0, mu_0, alpha_0) = do_analytic_backprop(
        sigma_mid, APL_g, ref.beta_gamma, ref.rigidity_Tm
    )

    # Longitudinal parameters (sigpt_0 [-] from input arguments)
    sigt_0 = 1e-3  # [m]
    emit_t = math.sqrt(sigt_0**2 * sigpt_0**2 - 0**2)
    print(
        f"sigt_0 = {sigt_0} [m], sigpt_0 = {sigpt_0} [-], emit_t = {emit_t}", flush=True
    )
    print(flush=True)

    #   particle bunch
    distr = distribution.Gaussian(
        lambdaX=math.sqrt(emitg / gamma_0),
        lambdaY=math.sqrt(emitg / gamma_0),
        lambdaT=sigt_0,  # OK for mutpt=0
        lambdaPx=math.sqrt(emitg / beta_0),
        lambdaPy=math.sqrt(emitg / beta_0),
        lambdaPt=sigpt_0,  # OK for mutpt=0
        muxpx=mu_0,
        muypy=mu_0,
        mutpt=0.0,
    )
    npart = 100000  # number of macro particles
    sim.add_particles(bunch_charge_C, distr, npart)

    # create the accelerator lattice

    # Plasma lens parameters for ConstF
    APL_k = APL_g / ref.rigidity_Tm
    APL_k_sqrt = np.sign(APL_k) * np.sqrt(np.abs(APL_k))
    print(f"APL_g = {APL_g} [T/m], APL_k = {APL_k} [1/m^2]", flush=True)
    print(flush=True)

    ns = 40  # number of slices per ds in the element
    monitor = elements.BeamMonitor("monitor", backend="h5")
    APL = None
    if lensType == "ChrPlasmaLens":
        APL = elements.ChrPlasmaLens(
            name="APL", ds=APL_length, k=APL_g, unit=1, nslice=ns
        )

    elif lensType == "ConstF":
        APL = elements.ConstF(
            name="APL", ds=APL_length, kx=APL_k_sqrt, ky=APL_k_sqrt, kt=0.0, nslice=ns
        )

    elif lensType == "ChrDrift":
        # For comparison with k=0
        assert float(APL_g) == 0.0
        APL = elements.ChrDrift(name="APL", ds=APL_length, nslice=ns)

    elif lensType == "ChrQuad":
        # For comparison with k != 0, single plane
        APL = elements.ChrQuad(name="APL", ds=APL_length, k=APL_g, unit=1, nslice=ns)
    else:
        raise ValueError(f"Unknown lensType {lensType}")

    lattice = [
        monitor,
        APL,
        monitor,
    ]

    # run simulation
    sim.lattice.extend(lattice)
    sim.track_particles()

    # clean shutdown
    sim.finalize()


def run_APL_envelope(
    APL_g: float, sigpt_0: float, sigma_mid: float, lensType: str = "ChrPlasmaLens"
):
    """
    Run a plasma lens envelope simulation with the given APL gradient APL_g [T/m], sigma_pt [-], and sigma_mid [m].
    Can specify lensType='ChrPlasmaLens' | 'ConstF' | 'ChrDrift' (expects APL_g = 0) | 'ChrQuad' (only horizontal plane valid)
    Note that for lensType='ChrPlasmaLens' | 'ChrDrift' | 'ChrQuad', the envelope simulation is not yet implemented.
    """

    from impactx import ImpactX, distribution, elements

    print("", flush=True)
    print(
        f"*** run_APL_envelope({APL_g}, {sigpt_0}, {sigma_mid}, '{lensType}') :",
        flush=True,
    )
    print("", flush=True)

    sim = ImpactX()

    # set numerical parameters and IO control
    sim.space_charge = False
    # sim.diagnostics = False  # benchmarking
    sim.slice_step_diagnostics = True

    # domain decomposition & space charge mesh
    sim.init_grids()

    ## Physics parameters for test (APL_g from input arguments)

    # Load a 200 MeV electron beam with alpha=0 (x and y)
    #  in the center of the APL
    # reference particle
    ref = sim.particle_container().ref_particle()
    ref.set_charge_qe(-1.0).set_mass_MeV(mass_MeV).set_kin_energy_MeV(kin_energy_MeV)

    (emitg, beta_0, gamma_0, mu_0, alpha_0) = do_analytic_backprop(
        sigma_mid, APL_g, ref.beta_gamma, ref.rigidity_Tm
    )

    # Longitudinal parameters (sigpt_0 [-] from input arguments)
    sigt_0 = 1e-3  # [m]
    emit_t = math.sqrt(sigt_0**2 * sigpt_0**2 - 0**2)
    print(
        f"sigt_0 = {sigt_0} [m], sigpt_0 = {sigpt_0} [-], emit_t = {emit_t}", flush=True
    )
    print(flush=True)

    #   particle bunch
    distr = distribution.Gaussian(
        lambdaX=math.sqrt(emitg / gamma_0),
        lambdaY=math.sqrt(emitg / gamma_0),
        lambdaT=sigt_0,  # OK for mutpt=0
        lambdaPx=math.sqrt(emitg / beta_0),
        lambdaPy=math.sqrt(emitg / beta_0),
        lambdaPt=sigpt_0,  # OK for mutpt=0
        muxpx=mu_0,
        muypy=mu_0,
        mutpt=0.0,
    )
    # npart = 100000  # number of macro particles
    # sim.add_particles(bunch_charge_C, distr, npart)

    # create the accelerator lattice

    # Plasma lens parameters for ConstF
    APL_k = APL_g / ref.rigidity_Tm
    APL_k_sqrt = np.sign(APL_k) * np.sqrt(np.abs(APL_k))
    print(f"APL_g = {APL_g} [T/m], APL_k = {APL_k} [1/m^2]", flush=True)
    print(flush=True)

    ns = 40  # number of slices per ds in the element
    monitor = elements.BeamMonitor("monitor", backend="h5")
    APL = None
    if lensType == "ChrPlasmaLens":
        APL = elements.ChrPlasmaLens(
            name="APL", ds=APL_length, k=APL_g, unit=1, nslice=ns
        )

    elif lensType == "ConstF":
        APL = elements.ConstF(
            name="APL", ds=APL_length, kx=APL_k_sqrt, ky=APL_k_sqrt, kt=0.0, nslice=ns
        )

    elif lensType == "ChrDrift":
        # For comparison with k=0
        assert float(APL_g) == 0.0
        APL = elements.ChrDrift(name="APL", ds=APL_length, nslice=ns)

    elif lensType == "ChrQuad":
        # For comparison with k != 0, single plane
        APL = elements.ChrQuad(name="APL", ds=APL_length, k=APL_g, unit=1, nslice=ns)
    else:
        raise ValueError(f"Unknown lensType {lensType}")

    lattice = [
        monitor,
        APL,
        monitor,
    ]

    # run simulation
    sim.lattice.extend(lattice)
    sim.init_envelope(ref, distr)
    sim.track_envelope()

    # clean shutdown
    sim.finalize()


def get_beta_gamma(Ek: float, m0: float):
    "Get beta*gamma given Ek and m_0c^2 in the same units (e.g. MeV)"
    E0 = Ek + m0
    gamma_rel = E0 / m0
    beta_rel = np.sqrt(
        (1.0 - 1.0 / np.sqrt(gamma_rel)) * (1.0 + 1.0 / np.sqrt(gamma_rel))
    )
    return beta_rel * gamma_rel


def get_rigidity(P0: float):
    "Convert momentum P0 [eV/c] to magnetic rigidity [T*m] for electron"
    return -P0 / 299792458


def do_analytic_backprop(sigma_mid: float, APL_g: float, beta_gamma: float, rigidity):
    """
    Make analytical back-propagation from lens midpoint to entry,
    taking the beam sigma at the midpoint [m] and APL gradient [T/m],
    along with the assumed beta_gamma and rigitdity [T*m].

    Also print analytical estimates of initial, midpoint, and final Twiss parameters;
    these can be used to compare simulation parameters with for tests.

    Returns:
      (emitg, beta_0, gamma_0, mu_0, alpha_0),
      which defines the beam parameters at the lens entry.
    """

    # Midpoint parameters
    alpha_mid = 0.0
    # sigma_mid = 10e-6  # [m]
    emitn = 10e-6  # [m]
    emitg = emitn / beta_gamma  # [m]
    beta_mid = sigma_mid**2 / emitg  # [m]
    gamma_mid = 1 / beta_mid  # [1/m]
    print(
        f"sigma_mid = {sigma_mid} [m], beta_mid = {beta_mid} [m], gamma_mid = {gamma_mid} [m], alpha_mid = {alpha_mid}",
        flush=True,
    )
    print(
        f"emitn = {emitn} [m], emitg = {emitg} [m], beta_gamma = {beta_gamma}, rigidity = {rigidity} [T*m]",
        flush=True,
    )
    print(flush=True)

    # Back-propagate 1/2 lens length as in vacuum,
    # from symmetry point in the middle of the lens to the start of the lens
    assert alpha_mid == 0.0
    beta_0 = beta_mid + (APL_length / 2) ** 2 / beta_mid
    alpha_0 = +APL_length / 2 / beta_mid
    gamma_0 = gamma_mid
    sigma_0 = math.sqrt(emitg * beta_0)
    sigmap_0 = math.sqrt(emitg * gamma_0)
    mu_0 = alpha_0 / math.sqrt(beta_0 * gamma_0)
    print(
        f"sigma_0 = {sigma_0} [m], beta_0 = {beta_0} [m], alpha_0 = {alpha_0}, sigmap_0 = {sigmap_0}",
        flush=True,
    )
    print(flush=True)

    # Forward-propagate that through the focusing/defocusing lens
    # from the beginning, ignoring energy spread
    (beta_end, alpha_end, gamma_end) = analytic_final_estimate(
        APL_g, rigidity, APL_length, beta_0, alpha_0
    )

    print(
        f"beta_end = {beta_end} [m], alpha_end = {alpha_end} [-], gamma_end = {gamma_end} [1/m]",
        flush=True,
    )
    sigma_end = np.sqrt(emitg * beta_end)
    sigmap_end = np.sqrt(emitg * gamma_end)
    print(f"sigma_end = {sigma_end} [m], sigmap_end = {sigmap_end} [-]", flush=True)
    print(flush=True)

    return (emitg, beta_0, gamma_0, mu_0, alpha_0)


def analytic_final_estimate(
    APL_g: float,
    rigidity: float,
    APL_length: float,
    beta_0: float,
    alpha_0: float,
    printk=True,
):
    "Analytical estimates of the beam Twiss parameters after the Plasma Lens"
    k = APL_g / rigidity
    if printk:
        print(f"k = {k} [1/m^2]", flush=True)
        print("", flush=True)

    if k > 0:
        M = np.asarray(
            [
                [
                    np.cos(APL_length * np.sqrt(k)),
                    np.sin(APL_length * np.sqrt(k)) / np.sqrt(k),
                ],
                [
                    -np.sqrt(k) * np.sin(APL_length * np.sqrt(k)),
                    np.cos(APL_length * np.sqrt(k)),
                ],
            ]
        )
    elif k < 0:
        M = np.asarray(
            [
                [
                    np.cosh(APL_length * np.sqrt(-k)),
                    np.sinh(APL_length * np.sqrt(-k)) / np.sqrt(-k),
                ],
                [
                    np.sqrt(-k) * np.sinh(APL_length * np.sqrt(-k)),
                    np.cosh(APL_length * np.sqrt(-k)),
                ],
            ]
        )
    else:
        M = np.asarray([[1, APL_length], [0, 1]])
    # Do the Twiss propagation
    B0 = np.asarray([[beta_0, -alpha_0], [-alpha_0, (1 + alpha_0**2) / beta_0]])
    B = M @ B0 @ M.T
    # print(B, flush=True)

    beta_end = B[0, 0]
    alpha_end = -B[0, 1]
    gamma_end = B[1, 1]

    return (beta_end, alpha_end, gamma_end)


def analytic_sigma_function(APL_g: float, sigma_mid: float):
    """
    Do a plasma lens analytical envelope calculation with
    the given APL gradient APL_g [T/m] and sigma_mid [m].
    """

    E_MeV = kin_energy_MeV + mass_MeV
    P0_MeV_c = np.sqrt(E_MeV**2 - mass_MeV**2)

    beta_gamma = get_beta_gamma(E_MeV, mass_MeV)
    rigidity = get_rigidity(P0_MeV_c * 1e6)  # [T*m]

    (emitg, beta_0, gamma_0, mu_0, alpha_0) = do_analytic_backprop(
        sigma_mid, APL_g, beta_gamma, rigidity
    )

    s = np.linspace(0, APL_length)
    sigma = np.empty_like(s)

    for i, s_ in enumerate(s):
        (beta_i, alpha_i, gamma_i) = analytic_final_estimate(
            APL_g, rigidity, s_, beta_0, alpha_0, printk=False
        )
        sigma[i] = np.sqrt(beta_i * emitg)

    print(s)
    print(sigma)

    return (s, sigma)
