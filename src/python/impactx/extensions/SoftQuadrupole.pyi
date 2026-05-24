"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Chad Mitchell, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

__all__: list[str] = ["register_SoftQuadrupole_extension"]

def register_SoftQuadrupole_extension(cls):
    """
    Extend SoftQuadrupole with an alternative constructor that accepts
    raw on-axis field/gradient data and computes Fourier coefficients
    internally.

    Parameters
    ----------
    cls : type
        The pybind11 ``SoftQuadrupole`` class to extend.
    """
