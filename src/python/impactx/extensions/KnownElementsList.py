"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""


def register_KnownElementsList_extension(kel):
    """KnownElementsList helper methods"""
    from ..madx_to_impactx import read_lattice
    from ..plot.Survey import plot_survey

    # register member functions for KnownElementsList
    kel.load_file = lambda self, madx_file, nslice=1: self.extend(
        read_lattice(madx_file, nslice)
    )
    kel.plot_survey = plot_survey
