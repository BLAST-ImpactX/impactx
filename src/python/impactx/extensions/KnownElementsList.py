"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
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
        from pals_schema.Line import Line

        # example: fodo.pals.yaml
        if extension == ".json":
            import json

            # Read the JSON data from the test file
            with open(filename, "r") as file:
                json_data = json.loads(file.read())
            # Parse the JSON data back into a Line object
            # (data validation happens here automatically)
            self.from_pals(Line(**json_data), nslice)
            return

        elif extension == ".yaml":
            import yaml

            # Read the YAML data from the test file
            with open(filename, "r") as file:
                yaml_data = yaml.safe_load(file)
            # Parse the YAML data back into a Line object
            # (data validation happens here automatically)
            self.from_pals(Line(**yaml_data), nslice)
            return

        # TODO: toml, xml

    raise RuntimeError(
        f"load_file: No support for file {filename} with extension {extension} yet."
    )


def from_pals(self, pals_line, nslice=1):
    """Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) Python Line.

    https://github.com/campa-consortium/pals-python
    """
    from pals_schema.DriftElement import DriftElement
    from pals_schema.QuadrupoleElement import QuadrupoleElement

    # Loop over the pals_line and create a new ImpactX KnownElementsList from it.
    #       Use self.extend(...) on the latter.
    print(nslice)
    ix_line = []
    for pals_element in pals_line.line:
        if isinstance(pals_element, DriftElement):
            ix_line.append(
                elements.Drift(
                    name=pals_element.name, ds=pals_element.length, nslice=nslice
                )
            )
        elif isinstance(pals_element, QuadrupoleElement):
            ix_line.append(
                elements.ChrQuad(
                    name=pals_element.name,
                    ds=pals_element.length,
                    k=pals_element.MagneticMultipoleP.Bn1,
                    unit=0,
                    nslice=nslice,
                )
            )

    self.extend(ix_line)


def register_KnownElementsList_extension(kel):
    """KnownElementsList helper methods"""
    from ..plot.Survey import plot_survey

    # register member functions for KnownElementsList
    kel.from_pals = from_pals
    kel.load_file = load_file
    kel.plot_survey = plot_survey
