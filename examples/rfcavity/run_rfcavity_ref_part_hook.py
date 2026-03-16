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


def get_effective_fourier_coeffs_AB(cos_coeffs, sin_coeffs, f, L, beta):
    w = 2.0 * np.pi * f
    k = w / sconst.c
    sinkL = np.sin(k * L / (2.0 * beta))

    Abasearr = np.array(
        [
            (-1) ** j
            / ((k * L - 2 * j * np.pi * beta) * (k * L + 2 * j * np.pi * beta))
            for j in range(ncoef)
        ]
    )
    Bbasearr = np.array(
        [
            (-1) ** j
            * j
            / ((k * L - 2 * j * np.pi * beta) * (k * L + 2 * j * np.pi * beta))
            for j in range(ncoef)
        ]
    )
    Acoeffs = -2.0 * Abasearr * k * L**2 * beta * cos_coeffs * sinkL
    Bcoeffs = 4.0 * Bbasearr * np.pi * L * beta**2 * sin_coeffs * sinkL
    Asum = sum(Acoeffs[1 : ncoef - 1]) + Acoeffs[0] / 2.0
    Bsum = sum(Bcoeffs[1 : ncoef - 1])

    A = Asum
    B = Bsum

    return A, B


def get_synch_phase_Veff(cos_coeffs, sin_coeffs, f, L, emax, beta, phase, t0):
    # predicted energy gain
    zmin = -L / (2.0 * beta)
    theta = phase * np.pi / 180.0
    w = 2.0 * np.pi * f
    k = w / sconst.c

    A, B = get_effective_fourier_coeffs_AB(cos_coeffs, sin_coeffs, f, L, beta)

    dpt = emax * (
        A * np.cos(k * (t0 - zmin / beta) + theta)
        + B * np.sin(k * (t0 - zmin / beta) + theta)
    )
    dgamma = -dpt

    Veff = emax * np.sqrt(A**2 + B**2)

    numerator = emax * (
        -B * np.cos(k * (t0 - zmin / beta) + theta)
        + A * np.sin(k * (t0 - zmin / beta) + theta)
    )

    synch_phase = np.arctan(numerator / dpt)

    return synch_phase, Veff, dgamma


def get_phase_emax(cos_coeffs, sin_coeffs, f, L, Veff, beta, synch_phase, t0):
    zmin = -L / (2.0 * beta)
    w = 2.0 * np.pi * f
    k = w / sconst.c

    A, B = get_effective_fourier_coeffs_AB(cos_coeffs, sin_coeffs, f, L, beta)

    term1 = -B - A * np.tan(synch_phase)
    term2 = -A + B * np.tan(synch_phase)
    phase = -k * (t0 - zmin / beta) + np.arctan2(term1, term2)
    phase_deg = np.degrees(np.mod(phase, 2 * np.pi))

    emax = Veff / np.sqrt(A**2 + B**2)

    return phase_deg, emax


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
        print(
            f"  Beam at s={ref.s:.2f}m, t={ref.t:.2f}s, gamma at entry={ref.gamma:.2f}",
            flush=True,
        )
        # Extract the effective RF voltage and synchronous phase from the ImpactX RFCavity input parameters:
        synch_phase, Veff, dgamma = get_synch_phase_Veff(
            cos_coeffs, sin_coeffs, f, L, emax, beta, phase, t0
        )
        print(
            f"  RF cavity synchronous phase={synch_phase:.2f}, V effective={Veff:.2f}, predicted change in gamma={dgamma:.2f}",
            flush=True,
        )
        # From the effective RF voltage and synchronous phase, set the ImpactX RFCavity input parameters:
        updated_phase, updated_emax = get_phase_emax(
            cos_coeffs, sin_coeffs, f, L, Veff, beta, synch_phase, t0
        )
        element.phase = updated_phase
        element.escale = updated_emax
        print(
            f"  RF cavity updated (reset) values of phase={updated_phase:.2f} and emax={updated_emax:.2f}",
            flush=True,
        )


sim.hook["before_element"] = hook_before_element

# run simulation
sim.track_reference(ref)

# clean shutdown
sim.finalize()
