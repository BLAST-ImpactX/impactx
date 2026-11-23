#!/usr/bin/env python3
#
# Copyright 2022-2026 ImpactX contributors
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""Shared vocabulary for the ImpactX element model tiers.

ImpactX implements most beamline elements at up to three levels of physical
fidelity, documented in ``docs/source/theory/assumptions.rst``:

============ =========================== ==========================
tier         class name prefix           example
============ =========================== ==========================
``linear``   (none)                      ``Drift``, ``Quad``
``paraxial`` ``Chr``                     ``ChrDrift``, ``ChrQuad``
``exact``    ``Exact``                   ``ExactDrift``, ``ExactQuad``
============ =========================== ==========================

Cheaper tiers run faster, richer tiers model more physics. This module holds
the one place of truth for that ordering, so that every user-facing knob that
selects a model speaks the same language. Those knobs are
``KnownElementsList.load_file(min_model=...)`` and
``FilteredElementsList.replace_with_drifts(model=...)``.

Not every element family implements every tier. There is no ``ChrSbend`` and no
``ExactSol``, for instance. :func:`select_model` therefore takes a table of the
tiers that a family *does* implement and rounds a requested minimum tier *up* to
the cheapest implemented one that satisfies it.
"""

from functools import partial

from impactx import elements

#: The model tiers, ordered from cheapest/simplest to most faithful.
MODEL_TIERS = ("linear", "paraxial", "exact")

_MODEL_RANK = {tier: rank for rank, tier in enumerate(MODEL_TIERS)}


def tier_rank(model: str) -> int:
    """Return the ordering rank of a model tier (higher = more faithful)."""
    return _MODEL_RANK[model]


def validate_model(model: str, *, argument: str = "min_model", extra_values=()) -> str:
    """Raise ``ValueError`` unless ``model`` names a known tier.

    :param model: the value to check
    :param argument: name of the user-facing argument, used in the error message
    :param extra_values: additional accepted values beyond :data:`MODEL_TIERS`
    :return: ``model``, unchanged
    """
    allowed = set(MODEL_TIERS) | set(extra_values)
    if model not in allowed:
        raise ValueError(f"{argument} must be one of {sorted(allowed)}, got {model!r}")
    return model


def tier_of_class(type_name: str) -> str:
    """Return the model tier of an ImpactX element class name.

    The class-name prefix *is* the tier. ``Exact*`` is exact, ``Chr*`` is
    paraxial, everything else is linear.
    """
    if type_name.startswith("Exact"):
        return "exact"
    if type_name.startswith("Chr"):
        return "paraxial"
    return "linear"


def select_model(builders: dict, min_model: str):
    """Pick the cheapest implemented model tier that is at least ``min_model``.

    :param builders: maps tier name to a callable creating the element. A family
        registers only the tiers ImpactX actually implements. A feature-driven
        requirement, such as a thick octupole needing the exact model, is
        expressed by omitting the cheaper tiers.
    :param min_model: the requested fidelity floor, one of :data:`MODEL_TIERS`
    :return: ``(tier, builder)``. If no implemented tier reaches ``min_model``,
        the highest implemented tier is returned instead. Callers detect that
        case by comparing the returned tier against ``min_model`` and warn.
    """
    floor = _MODEL_RANK[min_model]
    tiers = sorted(builders, key=_MODEL_RANK.__getitem__)
    for tier in tiers:
        if _MODEL_RANK[tier] >= floor:
            return tier, builders[tier]
    return tiers[-1], builders[tiers[-1]]


#: Drift models per tier. All three share one constructor signature.
DRIFT_MODEL_CLASSES = {
    "linear": elements.Drift,
    "paraxial": elements.ChrDrift,
    "exact": elements.ExactDrift,
}

#: Quadrupole models per tier. ``unit=0`` is the MAD-X convention (``k`` in
#: m^(-2)). Plain ``Quad`` is the only one without a ``unit`` argument and
#: always assumes that convention.
QUAD_MODEL_CLASSES = {
    "linear": elements.Quad,
    "paraxial": partial(elements.ChrQuad, unit=0),
    "exact": partial(elements.ExactQuad, unit=0),
}
