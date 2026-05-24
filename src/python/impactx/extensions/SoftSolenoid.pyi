"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Chad Mitchell, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

__all__: list[str] = ["register_SoftSolenoid_extension"]

def register_SoftSolenoid_extension(cls):
    """
    Extend SoftSolenoid with an alternative constructor that accepts
    raw on-axis field data and computes Fourier coefficients internally.

    Parameters
    ----------
    cls : type
        The pybind11 ``SoftSolenoid`` class to extend.
    """
