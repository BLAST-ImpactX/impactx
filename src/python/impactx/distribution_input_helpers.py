#!/usr/bin/env python3
#
# Copyright 2024 ImpactX contributors
# Authors: Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np


def twiss(
    beta_x: np.float64,
    beta_y: np.float64,
    beta_t: np.float64,
    emitt_x: np.float64,
    emitt_y: np.float64,
    emitt_t: np.float64,
    alpha_x: np.float64 = 0.0,
    alpha_y: np.float64 = 0.0,
    alpha_t: np.float64 = 0.0,
    meanX: np.float64 = 0.0,
    meanY: np.float64 = 0.0,
    meanT: np.float64 = 0.0,
    meanPx: np.float64 = 0.0,
    meanPy: np.float64 = 0.0,
    meanPt: np.float64 = 0.0,
    dispX: np.float64 = 0.0,
    dispY: np.float64 = 0.0,
    dispPx: np.float64 = 0.0,
    dispPy: np.float64 = 0.0,
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
    :param meanX: offset of the mean (centroid) position in x from that of the reference particle
    :param meanY: offset of the mean (centroid) position in y from that of the reference particle
    :param meanT: offset of the mean (centroid) position in t from that of the reference particle
    :param meanPx: offset of the mean (centroid) momentum in x from that of the reference particle
    :param meanPy: offset of the mean (centroid) momentum in y from that of the reference particle
    :param meanPt: offset of the mean (centroid) momentum in t from that of the reference particle
    :param dispX: dispersion and its derivative in horizontal and vertical directions
    :param dispY: dispersion and its derivative in horizontal and vertical directions
    :param dispPx: dispersion and its derivative in horizontal and vertical directions
    :param dispPy: dispersion and its derivative in horizontal and vertical directions
    :return: A dictionary containing calculated phase space input: 'lambdaX', 'lambdaY', 'lambdaT', 'lambdaPx', 'lambdaPy', 'lambdaPt', 'muxpx', 'muypy', 'mutpt'.
    """
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
        "meanX": meanX,
        "meanY": meanY,
        "meanT": meanT,
        "meanPx": meanPx,
        "meanPy": meanPy,
        "meanPt": meanPt,
        "dispX": dispX,
        "dispY": dispY,
        "dispPx": dispPx,
        "dispPy": dispPy,
    }
