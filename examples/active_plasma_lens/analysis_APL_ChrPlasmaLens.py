#!/usr/bin/env python3
#
# Copyright 2022-2025 ImpactX contributors
# Authors: Axel Huebl, Chad Mitchell, Kyrre Sjobak
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import openpmd_api as io
from scipy.stats import moment


def get_moments(beam):
    """Calculate standard deviations of beam position & momenta
    and emittance values

    Returns
    -------
    sigx, sigy, sigt, emittance_x, emittance_y, emittance_t
    """
    sigx = moment(beam["position_x"], moment=2) ** 0.5  # variance -> std dev.
    sigpx = moment(beam["momentum_x"], moment=2) ** 0.5
    sigy = moment(beam["position_y"], moment=2) ** 0.5
    sigpy = moment(beam["momentum_y"], moment=2) ** 0.5
    sigt = moment(beam["position_t"], moment=2) ** 0.5
    sigpt = moment(beam["momentum_t"], moment=2) ** 0.5

    epstrms = beam.cov(ddof=0)
    emittance_x = (sigx**2 * sigpx**2 - epstrms["position_x"]["momentum_x"] ** 2) ** 0.5
    emittance_y = (sigy**2 * sigpy**2 - epstrms["position_y"]["momentum_y"] ** 2) ** 0.5
    emittance_t = (sigt**2 * sigpt**2 - epstrms["position_t"]["momentum_t"] ** 2) ** 0.5

    return (sigx, sigy, sigt, emittance_x, emittance_y, emittance_t)


def get_twiss(beam):
    "Calculate the beam Twiss parameters from position and momenta values"

    epstrms = beam.cov(ddof=0)

    sigx2 = epstrms["position_x"]["position_x"]
    sigpx2 = epstrms["momentum_x"]["momentum_x"]
    emittance_x = (sigx2 * sigpx2 - epstrms["position_x"]["momentum_x"] ** 2) ** 0.5
    beta_x = sigx2 / emittance_x
    alpha_x = -epstrms["position_x"]["momentum_x"] / emittance_x

    sigy2 = epstrms["position_y"]["position_y"]
    sigpy2 = epstrms["momentum_y"]["momentum_y"]
    emittance_y = (sigy2 * sigpy2 - epstrms["position_y"]["momentum_y"] ** 2) ** 0.5
    beta_y = sigy2 / emittance_y
    alpha_y = -epstrms["position_y"]["momentum_y"] / emittance_y

    return (beta_x, beta_y, alpha_x, alpha_y)


def get_beams():
    "Load the initial and final beam from last simulation"
    series = io.Series("diags/openPMD/monitor.h5", io.Access.read_only)
    last_step = list(series.iterations)[-1]
    initial = series.iterations[1].particles["beam"].to_df()
    beam_final = series.iterations[last_step].particles["beam"]
    final = beam_final.to_df()

    return (initial, beam_final, final)
