"""
This file is part of ImpactX

Copyright 2022-2025 ImpactX contributors
Authors: Chad Mitchell, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

import numpy as np

__all__: list[str] = ["fourier_coefficients", "np"]

def fourier_coefficients(z, field_or_gradient, ncoef):
    """
    Calculate Fourier coefficients of on-axis field data.

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
