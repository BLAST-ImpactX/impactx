"""
Compensated arithmetic utilities for ImpactX.

This module provides high-precision arithmetic functions using the CompensatedParticleReal
type, which implements the Klein (2006) second-order iterative Kahan-Babuška algorithm.

Copyright 2022-2024 The ImpactX Community
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""


def compensated_sum(values):
    """
    Compute the sum of values using compensated arithmetic.

    This function uses CompensatedParticleReal to maintain high precision
    when summing many small values, avoiding floating-point precision errors.

    Parameters
    ----------
    values : iterable
        An iterable of numeric values to sum

    Returns
    -------
    float
        The compensated sum of all values

    Examples
    --------
    >>> compensated_sum([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    1.0
    """
    from .impactx_pybind import CompensatedParticleReal

    result = CompensatedParticleReal(0.0)
    for value in values:
        result += value
    return result.value


def compensated_cumsum(values):
    """
    Compute the cumulative sum of values using compensated arithmetic.

    This function uses CompensatedParticleReal to maintain high precision
    when computing cumulative sums of many small values, avoiding floating-point
    precision errors that can accumulate in long sequences.

    Parameters
    ----------
    values : iterable
        An iterable of numeric values to compute cumulative sum for

    Returns
    -------
    list
        A list containing the cumulative sum at each position, starting with 0.0

    Examples
    --------
    >>> compensated_cumsum([0.1, 0.2, 0.3])
    [0.0, 0.1, 0.3, 0.6]
    """
    from .impactx_pybind import CompensatedParticleReal

    result = CompensatedParticleReal(0.0)
    cumulative = [0.0]  # Start with initial value
    for value in values:
        result += value
        cumulative.append(result.value)
    return cumulative
