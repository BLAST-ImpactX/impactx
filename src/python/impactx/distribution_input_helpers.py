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
    mean_x: np.float64 = 0.0,
    mean_y: np.float64 = 0.0,
    mean_t: np.float64 = 0.0,
    mean_px: np.float64 = 0.0,
    mean_py: np.float64 = 0.0,
    mean_pt: np.float64 = 0.0,
    dispersion_x: np.float64 = 0.0,
    dispersion_y: np.float64 = 0.0,
    dispersion_px: np.float64 = 0.0,
    dispersion_py: np.float64 = 0.0,
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
    :param dispersion_x: dispersion and its derivative in horizontal and vertical directions
    :param dispersion_y: dispersion and its derivative in horizontal and vertical directions
    :param dispersion_px: dispersion and its derivative in horizontal and vertical directions
    :param dispersion_py: dispersion and its derivative in horizontal and vertical directions
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
        "lambda_x": np.sqrt(emitt_x / gamma_x),
        "lambda_y": np.sqrt(emitt_y / gamma_y),
        "lambda_t": np.sqrt(emitt_t / gamma_t),
        "lambda_px": np.sqrt(emitt_x / beta_x),
        "lambda_py": np.sqrt(emitt_y / beta_y),
        "lambda_pt": np.sqrt(emitt_t / beta_t),
        "mu_x_px": alpha_x / np.sqrt(beta_x * gamma_x),
        "mu_y_py": alpha_y / np.sqrt(beta_y * gamma_y),
        "mu_t_pt": alpha_t / np.sqrt(beta_t * gamma_t),
        "mean_x": mean_x,
        "mean_y": mean_y,
        "mean_t": mean_t,
        "mean_px": mean_px,
        "mean_py": mean_py,
        "mean_pt": mean_pt,
        "dispersion_x": dispersion_x,
        "dispersion_y": dispersion_y,
        "dispersion_px": dispersion_px,
        "dispersion_py": dispersion_py,
    }
