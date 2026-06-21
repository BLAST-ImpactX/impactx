"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell, Edoardo Zoni
License: BSD-3-Clause-LBNL
"""

from impactx import elements


def flatten_pals(pals_data, registry=None):
    """Flatten a PALS root, lattice, or beamline to a list of PALS elements.

    Placeholder references are resolved from the root facility definitions.
    """
    from pals import BeamLine, Lattice, PALSroot, PlaceholderName

    if registry is None:
        registry = {}

    if isinstance(pals_data, PALSroot):
        registry = {
            item.name: item
            for item in pals_data.facility
            if not isinstance(item, PlaceholderName) and hasattr(item, "name")
        }

        if len(pals_data.facility) == 1:
            return flatten_pals(pals_data.facility[0], registry)

        active_entry = pals_data.facility[-1]
        if not isinstance(active_entry, PlaceholderName):
            raise RuntimeError(
                "from_pals: PALS roots with multiple facility entries must "
                "select the active lattice or beamline with a final 'use' entry."
            )
        return flatten_pals(active_entry, registry)

    if isinstance(pals_data, PlaceholderName):
        if pals_data.element is not None:
            return flatten_pals(pals_data.element, registry)
        if pals_data.name not in registry:
            raise RuntimeError(
                f"from_pals: Cannot resolve PALS element reference {pals_data.name!r}."
            )
        return flatten_pals(registry[pals_data.name], registry)

    if isinstance(pals_data, Lattice):
        if len(pals_data.branches) != 1:
            raise RuntimeError(
                "from_pals: ImpactX currently supports PALS lattices with exactly "
                f"one branch, but got {len(pals_data.branches)}."
            )
        return flatten_pals(pals_data.branches[0], registry)

    if isinstance(pals_data, BeamLine):
        pals_elements = []
        for element in pals_data.line:
            pals_elements.extend(flatten_pals(element, registry))
        return pals_elements

    return [pals_data]


def read_lattice(pals_beamline, nslice=1):
    """Translate a Particle Accelerator Lattice Standard (PALS) object into ImpactX elements.

    https://github.com/campa-consortium/pals-python

    :param pals_beamline: a PALS root, lattice, or beamline object
    :param nslice: number of ds slices per element
    :return: list of ImpactX.KnownElements
    """
    from pals import Drift, Quadrupole

    pals_elements = flatten_pals(pals_beamline)

    # Loop over the flattened PALS beamline and create the matching ImpactX elements.
    ix_beamline = []
    for pals_element in pals_elements:
        if isinstance(pals_element, Drift):
            ix_beamline.append(
                elements.Drift(
                    name=pals_element.name, ds=pals_element.length, nslice=nslice
                )
            )
        elif isinstance(pals_element, Quadrupole):
            magnetic_multipole = pals_element.MagneticMultipoleP
            if magnetic_multipole is None:
                raise RuntimeError(
                    f"from_pals: No magnetic multipole input provided for element of kind {type(pals_element)}."
                )

            if getattr(magnetic_multipole, "Bn1", None) is not None:
                k_quad = magnetic_multipole.Bn1
                unit_quad = 1
            elif getattr(magnetic_multipole, "Kn1", None) is not None:
                k_quad = magnetic_multipole.Kn1
                unit_quad = 0
            else:
                raise RuntimeError(
                    f"from_pals: No gradient input provided for element of kind {type(pals_element)}."
                )
            ix_beamline.append(
                elements.ChrQuad(
                    name=pals_element.name,
                    ds=pals_element.length,
                    k=k_quad,
                    unit=unit_quad,
                    nslice=nslice,
                )
            )
        else:
            raise RuntimeError(
                f"from_pals: No support for elements of kind {type(pals_element)} yet."
            )

    return ix_beamline
