"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell, Edoardo Zoni
License: BSD-3-Clause-LBNL
"""

import os

from impactx import elements


def load_file(self, filename, nslice=1):
    """Load and append a lattice file from MAD-X (.madx) or PALS (e.g., .pals.yaml) formats."""

    # Attempt to strip two levels of file extensions to determine the schema.
    #   Examples: fodo.madx, fodo.pals.yaml, fodo.pals.json, ...
    file_noext, extension = os.path.splitext(filename)
    file_noext_noext, extension_inner = os.path.splitext(file_noext)

    if extension == ".madx":
        # example: fodo.madx
        from ..madx_to_impactx import read_lattice

        self.extend(read_lattice(filename, nslice))
        return

    elif extension_inner == ".pals":
        from pals_schema.BeamLine import BeamLine

        # examples: fodo.pals.yaml, fodo.pals.json
        with open(filename, "r") as file:
            if extension == ".json":
                import json

                pals_data = json.loads(file.read())
            elif extension == ".yaml":
                import yaml

                pals_data = yaml.safe_load(file)
            # TODO: toml, xml
            else:
                raise RuntimeError(
                    f"load_file: No support for PALS file {filename} with extension {extension} yet."
                )

        # Parse the data dictionary back into a PALS `BeamLine` object.
        # The automatically PALS data validation happens here.
        self.from_pals(BeamLine(**pals_data), nslice)
        return

    raise RuntimeError(
        f"load_file: No support for file {filename} with extension {extension} yet."
    )


def from_pals(self, pals_beamline, nslice=1):
    """Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) Python BeamLine.

    https://github.com/campa-consortium/pals-python
    """
    from pals_schema.DriftElement import DriftElement
    from pals_schema.QuadrupoleElement import QuadrupoleElement

    # Loop over the pals_beamline and create a new ImpactX KnownElementsList from it.
    #       Use self.extend(...) on the latter.
    ix_beamline = []
    for pals_element in pals_beamline.line:
        if isinstance(pals_element, DriftElement):
            ix_beamline.append(
                elements.Drift(
                    name=pals_element.name, ds=pals_element.length, nslice=nslice
                )
            )
        elif isinstance(pals_element, QuadrupoleElement):
            ix_beamline.append(
                elements.ChrQuad(
                    name=pals_element.name,
                    ds=pals_element.length,
                    k=pals_element.MagneticMultipoleP.Bn1,
                    unit=0,
                    nslice=nslice,
                )
            )
        else:
            raise RuntimeError(
                f"from_pals: No support for elements of kind {type(pals_element)} yet."
            )

    self.extend(ix_beamline)


def register_KnownElementsList_extension(kel):
    """KnownElementsList helper methods"""
    from ..plot.Survey import plot_survey

    # register member functions for KnownElementsList
    kel.from_pals = from_pals
    kel.load_file = load_file
    kel.plot_survey = plot_survey
