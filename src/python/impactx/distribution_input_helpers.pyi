from __future__ import annotations

import numpy as numpy

__all__: list[str] = ["numpy", "twiss"]

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
