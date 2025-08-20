"""

This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell, Edoardo Zoni
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

import os as os

from impactx.impactx_pybind import elements

__all__: list[str] = [
    "elements",
    "from_pals",
    "load_file",
    "os",
    "register_KnownElementsList_extension",
]

def from_pals(self, pals_beamline, nslice=1):
    """
    Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) Python BeamLine.

        https://github.com/campa-consortium/pals-python

    """

def load_file(self, filename, nslice=1):
    """
    Load and append a lattice file from MAD-X (.madx) or PALS (e.g., .pals.yaml) formats.
    """

def register_KnownElementsList_extension(kel):
    """
    KnownElementsList helper methods
    """
