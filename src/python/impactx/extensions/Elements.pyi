"""
This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl
License: BSD-3-Clause-LBNL

Adds value-based ``__eq__``, ``__hash__``, and ``isclose()`` to every
lattice element class that exposes a ``to_dict()`` method.

The semantics are derived from the dict returned by ``to_dict()``:

- ``__eq__`` compares dict keys and values via ``==`` (with element-wise
  comparison for lists and ``numpy.array_equal`` for AMReX SmallMatrix
  values). It returns ``NotImplemented`` for foreign operands so Python's
  reflected-equality fallback applies.
- ``__hash__`` is consistent with ``__eq__``. Hashing reflects the
  element's *current* parameter values; mutating an element after using
  it as a ``set`` member or ``dict`` key is undefined behavior, the same
  contract as a hypothetical hashable list.
- ``isclose(other, *, rtol, atol)`` mirrors ``math.isclose`` /
  ``numpy.isclose`` naming. Floats use ``math.isclose``; lists of floats
  are compared element-wise; AMReX SmallMatrix values use
  ``numpy.allclose``; everything else (ints, strings, ``None``) uses
  ``==``.
"""

from __future__ import annotations

import inspect as inspect
import math as math

__all__: list[str] = ["inspect", "isclose", "math", "register_elements_value_semantics"]

def _canon(v):
    """
    Map a ``to_dict()`` value to a hashable, type-tagged form.

    The type tag prevents accidental collisions between different value
    shapes (e.g., a list and a tuple with the same numeric contents).
    """

def _dict_close(d1, d2, rtol, atol): ...
def _dict_eq(d1, d2): ...
def _element_eq(self, other):
    """
    Value-based equality.

    Two elements are equal iff they are instances of the same class and
    their ``to_dict()`` outputs match key-for-key. Float values use
    ``==`` (so ``NaN != NaN``); list values are compared element-wise;
    AMReX SmallMatrix values use ``numpy.array_equal``.

    Returns ``NotImplemented`` when ``other`` is not the same element
    type, so Python's reflected-equality fallback applies (e.g., a
    foreign type's ``__eq__`` gets a chance, ultimately falling back to
    identity which yields ``False``).
    """

def _element_hash(self):
    """
    Value-based hash, consistent with ``__eq__``.

    Two elements that compare equal under ``==`` produce the same hash.
    The hash reflects the element's *current* parameter values; mutating
    an element after using it as a ``set`` member or ``dict`` key
    invalidates the container, the same contract a hashable mutable
    object would have.
    """

def _element_isclose(self, other, *, rtol=1e-12, atol=0.0, ignore_attributes=None):
    """
    Tolerant equality for lattice elements.

    Mirrors ``math.isclose`` / ``numpy.isclose`` naming. Float-valued
    fields are compared via ``math.isclose(rel_tol=rtol, abs_tol=atol)``;
    lists of floats are compared element-wise; AMReX SmallMatrix values
    use ``numpy.allclose``. All other value types (ints, strings,
    ``None``) fall back to strict ``==``. Mismatched element types and
    foreign operands return ``False`` (unless ``"type"`` is listed in
    ``ignore_attributes`` — see below).

    Parameters
    ----------
    other
        Element to compare against.
    rtol : float, optional
        Relative tolerance forwarded to ``math.isclose`` /
        ``numpy.allclose``. Default is ``1e-12`` — matches the
        ``dicts_equal`` helper used in the serialization tests, stricter
        than the ``math.isclose`` and ``numpy.isclose`` defaults.
    atol : float, optional
        Absolute tolerance. Default is ``0.0``.
    ignore_attributes : str or iterable of str, optional
        ``to_dict()`` keys to skip when comparing. Useful when comparing
        loaded files where some bookkeeping fields (``"name"``) or even
        the element variant (``"type"``) should not affect the verdict.

        Including ``"type"`` disables the strict same-class check, so
        e.g. ``Drift`` and ``ExactDrift`` can be compared on their
        common parameters. Any remaining keys must still match between
        the two ``to_dict()`` outputs.
    """

def _normalize_ignore(ignore_attributes):
    """
    Coerce ``ignore_attributes`` to a frozenset of key names.

    Accepts ``None``, a single string, or any iterable of strings.
    """

def _scalar_close(a, b, rtol, atol): ...
def _values_close(a, b, rtol, atol): ...
def _values_equal(a, b): ...
def isclose(a, b, **kwargs):
    """
    Free-function form of ``a.isclose(b)``.

    Provided so the call site reads symmetrically in ``a`` and ``b``,
    like :py:func:`math.isclose` / :py:func:`numpy.isclose`. Equivalent
    to the method form; forwards ``rtol``, ``atol``, and
    ``ignore_attributes``. Dispatches to ``b`` if ``a`` lacks
    ``isclose`` (e.g., a plain Python ``list`` paired with a
    ``KnownElementsList``).
    """

def register_elements_value_semantics(elements_module):
    """
    Attach ``__eq__``, ``__hash__``, and ``isclose`` to every class
    in ``elements_module`` that exposes a ``to_dict`` method, and the
    free-function ``isclose`` on the module itself.

    Parameters
    ----------
    elements_module : module
        Typically ``impactx.impactx_pybind.elements``.
    """

_DEFAULT_ATOL: float = 0.0
_DEFAULT_RTOL: float = 1e-12
_SKIP: frozenset
