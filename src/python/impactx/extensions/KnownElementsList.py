"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL
"""

import os


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
        import pals

        # example: fodo.pals.yaml
        if extension == ".json":
            import json

            # Read the JSON data from the test file
            with open(filename, "r") as file:
                json_data = json.loads(file.read())
            # Parse the JSON data back into a Line object
            self.from_pals(pals.Line(**json_data), nslice)
            return

        elif extension == ".yaml":
            import yaml

            # Read the JSON data from the test file
            with open(filename, "r") as file:
                yaml_data = yaml.safe_load(file)
            # Parse the JSON data back into a Line object
            self.from_pals(pals.Line(**yaml_data), nslice)
            return

        # TODO: toml, xml

    raise RuntimeError(
        f"load_file: No support for file {filename} with extension {extension} yet."
    )


def from_pals(self, pals_line, nslice=1):
    """Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) Python Line.

    https://github.com/campa-consortium/pals-python
    """
    # TODO: Loop over the pals_line and create a new ImpactX KnownElementsList from it.
    #       Use self.extend(...) on the latter.
    pass


def register_KnownElementsList_extension(kel):
    """KnownElementsList helper methods"""
    from ..plot.Survey import plot_survey

    # register member functions for KnownElementsList
    kel.from_pals = from_pals
    kel.load_file = load_file
    kel.plot_survey = plot_survey
