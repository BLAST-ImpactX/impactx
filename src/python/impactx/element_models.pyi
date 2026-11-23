"""
Shared vocabulary for the ImpactX element model tiers.

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

from __future__ import annotations

from functools import partial

import impactx.impactx_pybind.elements
from impactx.impactx_pybind import elements

__all__: list[str] = [
    "DRIFT_MODEL_CLASSES",
    "MODEL_TIERS",
    "QUAD_MODEL_CLASSES",
    "elements",
    "partial",
    "select_model",
    "tier_of_class",
    "tier_rank",
    "validate_model",
]

def select_model(builders: dict, min_model: str):
    """
    Pick the cheapest implemented model tier that is at least ``min_model``.

    :param builders: maps tier name to a callable creating the element. A family
        registers only the tiers ImpactX actually implements. A feature-driven
        requirement, such as a thick octupole needing the exact model, is
        expressed by omitting the cheaper tiers.
    :param min_model: the requested fidelity floor, one of :data:`MODEL_TIERS`
    :return: ``(tier, builder)``. If no implemented tier reaches ``min_model``,
        the highest implemented tier is returned instead. Callers detect that
        case by comparing the returned tier against ``min_model`` and warn.
    """

def tier_of_class(type_name: str) -> str:
    """
    Return the model tier of an ImpactX element class name.

    The class-name prefix *is* the tier. ``Exact*`` is exact, ``Chr*`` is
    paraxial, everything else is linear.
    """

def tier_rank(model: str) -> int:
    """
    Return the ordering rank of a model tier (higher = more faithful).
    """

def validate_model(
    model: str, *, argument: str = "min_model", extra_values=tuple()
) -> str:
    """
    Raise ``ValueError`` unless ``model`` names a known tier.

    :param model: the value to check
    :param argument: name of the user-facing argument, used in the error message
    :param extra_values: additional accepted values beyond :data:`MODEL_TIERS`
    :return: ``model``, unchanged
    """

DRIFT_MODEL_CLASSES: dict = {
    "linear": impactx.impactx_pybind.elements.Drift,
    "paraxial": impactx.impactx_pybind.elements.ChrDrift,
    "exact": impactx.impactx_pybind.elements.ExactDrift,
}
MODEL_TIERS: tuple = ("linear", "paraxial", "exact")
QUAD_MODEL_CLASSES: dict
_MODEL_RANK: dict = {"linear": 0, "paraxial": 1, "exact": 2}
