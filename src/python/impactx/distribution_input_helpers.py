#!/usr/bin/env python3
#
# Copyright 2024 ImpactX contributors
# Authors: Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy


def twiss(
    beta_x: numpy.float64,
    beta_y: numpy.float64,
    beta_t: numpy.float64,
    emitt_x: numpy.float64,
    emitt_y: numpy.float64,
    emitt_t: numpy.float64,
    alpha_x: numpy.float64 = 0.0,
    alpha_y: numpy.float64 = 0.0,
    alpha_t: numpy.float64 = 0.0,
    mean_x: numpy.float64 = 0.0,
    mean_y: numpy.float64 = 0.0,
    mean_t: numpy.float64 = 0.0,
    mean_px: numpy.float64 = 0.0,
    mean_py: numpy.float64 = 0.0,
    mean_pt: numpy.float64 = 0.0,
    dispersion_x: numpy.float64 = 0.0,
    dispersion_y: numpy.float64 = 0.0,
    dispersion_px: numpy.float64 = 0.0,
    dispersion_py: numpy.float64 = 0.0,
):
    """
    Helper function to convert Courant-Snyder / Twiss input into phase space ellipse input.

    :param beta_x: Beta function value (unit: meter) in the x dimension, must be a non-zero positive value.
    :param beta_y: Beta function value (unit: meter) in the y dimension, must be a non-zero positive value.
    :param beta_t: Beta function value (unit: meter) in the t dimension (arrival time differences multiplied by light speed), must be a non-zero positive value.
    :param emitt_x: Emittance value (unit: meter times radian) in the x dimension, must be a non-zero positive value.
    :param emitt_y: Emittance value (unit: meter times radian) in the y dimension, must be a non-zero positive value.
    :param emitt_t: Emittance value (unit: meter times radian) in the t dimension (arrival time differences multiplied by light speed), must be a non-zero positive value.
    :param alpha_x: Alpha function value () in the x dimension, default is 0.0.
    :param alpha_y: Alpha function value in the y dimension, default is 0.0.
    :param alpha_t: Alpha function value in the t dimension, default is 0.0.
    :param mean_x: offset of the mean (centroid) position in x from that of the reference particle
    :param mean_y: offset of the mean (centroid) position in y from that of the reference particle
    :param mean_t: offset of the mean (centroid) position in t from that of the reference particle
    :param mean_px: offset of the mean (centroid) momentum in x from that of the reference particle
    :param mean_py: offset of the mean (centroid) momentum in y from that of the reference particle
    :param mean_pt: offset of the mean (centroid) momentum in t from that of the reference particle
    :param dispersion_x: horizontal dispersion [m]
    :param dispersion_y: vertical dispersion [m]
    :param dispersion_px: derivative of horizontal dispersion
    :param dispersion_py: derivative of vertical dispersion
    :return: A dictionary containing calculated phase space input: 'lambdaX', 'lambdaY', 'lambdaT', 'lambdaPx', 'lambdaPy', 'lambdaPt', 'muxpx', 'muypy', 'mutpt'.
    """
    import numpy as np

    if beta_x <= 0.0 or beta_y <= 0.0 or beta_t <= 0.0:
        raise ValueError(
            "Input Error: The beta function values need to be non-zero positive values in all dimensions."
        )

    if emitt_x <= 0.0 or emitt_y <= 0.0 or emitt_t <= 0.0:
        raise ValueError(
            "Input Error: Emittance values need to be non-zero positive values in all dimensions."
        )

    betas = [beta_x, beta_y, beta_t]
    alphas = [alpha_x, alpha_y, alpha_t]

    gammas = []
    # calculate Courant-Snyder gammas
    for i in range(3):
        gammas.append((1 + alphas[i] ** 2) / betas[i])
    gamma_x, gamma_y, gamma_t = gammas

    return {
        "lambdaX": np.sqrt(emitt_x / gamma_x),
        "lambdaY": np.sqrt(emitt_y / gamma_y),
        "lambdaT": np.sqrt(emitt_t / gamma_t),
        "lambdaPx": np.sqrt(emitt_x / beta_x),
        "lambdaPy": np.sqrt(emitt_y / beta_y),
        "lambdaPt": np.sqrt(emitt_t / beta_t),
        "muxpx": alpha_x / np.sqrt(beta_x * gamma_x),
        "muypy": alpha_y / np.sqrt(beta_y * gamma_y),
        "mutpt": alpha_t / np.sqrt(beta_t * gamma_t),
        "meanX": mean_x,
        "meanY": mean_y,
        "meanT": mean_t,
        "meanPx": mean_px,
        "meanPy": mean_py,
        "meanPt": mean_pt,
        "dispX": dispersion_x,
        "dispY": dispersion_y,
        "dispPx": dispersion_px,
        "dispPy": dispersion_py,
    }
