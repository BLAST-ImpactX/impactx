"""
This file is part of ImpactX

Copyright 2022-2025 ImpactX contributors
Authors: Chad Mitchell, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import numpy as np


def fourier_coefficients(z, field_or_gradient, ncoef):
    """Calculate Fourier coefficients of on-axis field data.

    Uses the trapezoidal rule with linear interpolation to compute
    cosine and sine Fourier coefficients of the field profile,
    centered about the midpoint of the data range.

    Parameters
    ----------
    z : numpy.ndarray
        Longitudinal positions in meters, shape (N,), covering the
        element from entry (``min(z)``) to exit (``max(z)``).
        The range is scaled to the element length ``ds``.
    field_or_gradient : numpy.ndarray
        On-axis field or field gradient values, shape (N,),
        typically normalized to a peak absolute value of 1.
        These values are multiplied by the element's scaling
        parameter (``gscale``, ``bscale``, or ``escale``).
    ncoef : int
        Number of Fourier coefficients to compute.

    Returns
    -------
    cos_coeffs : numpy.ndarray
        Cosine Fourier coefficients, length *ncoef*.
    sin_coeffs : numpy.ndarray
        Sine Fourier coefficients, length *ncoef*.
    """
    ndatareal = len(z)

    zlen = z[-1] - z[0]
    zmid = (z[-1] + z[0]) / 2
    zhalf = zlen / 2.0
    h = zlen / (ndatareal - 1)

    j = np.arange(ncoef)  # (ncoef,)

    # Endpoint correction (trapezoidal rule)
    zz0 = z[0] - zmid
    zz1 = z[-1] - zmid
    angle0 = j * 2 * np.pi * zz0 / zlen  # (ncoef,)
    angle1 = j * 2 * np.pi * zz1 / zlen  # (ncoef,)

    cos_coeffs = (-0.5 * field_or_gradient[0] * np.cos(angle0) * h) / zhalf
    sin_coeffs = (-0.5 * field_or_gradient[0] * np.sin(angle0) * h) / zhalf
    cos_coeffs -= (0.5 * field_or_gradient[-1] * np.cos(angle1) * h) / zhalf
    sin_coeffs -= (0.5 * field_or_gradient[-1] * np.sin(angle1) * h) / zhalf

    # Interior points: interpolate field onto uniform grid, then integrate
    zz_uniform = np.arange(ndatareal) * h + z[0]  # (ndatareal,)
    ez1 = np.interp(zz_uniform, z, field_or_gradient)  # (ndatareal,)

    zz_centered = zz_uniform - zmid  # (ndatareal,)
    # Outer product: angles[i, j] = j * 2 * pi * zz_centered[i] / zlen
    angles = np.outer(zz_centered, j * 2 * np.pi / zlen)  # (ndatareal, ncoef)

    cos_coeffs += np.sum(ez1[:, np.newaxis] * np.cos(angles) * h, axis=0) / zhalf
    sin_coeffs += np.sum(ez1[:, np.newaxis] * np.sin(angles) * h, axis=0) / zhalf

    return cos_coeffs, sin_coeffs
