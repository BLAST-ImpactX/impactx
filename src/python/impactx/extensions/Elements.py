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

import inspect
import math

_DEFAULT_RTOL = 1e-12
_DEFAULT_ATOL = 0.0

# Names in the elements module that should not receive these methods.
# KnownElementsList / FilteredElementsList are containers, not elements;
# they do not currently expose a singular ``to_dict``, but we list them
# explicitly so the intent is unambiguous.
_SKIP = frozenset({"KnownElementsList", "FilteredElementsList"})


def _canon(v):
    """Map a ``to_dict()`` value to a hashable, type-tagged form.

    The type tag prevents accidental collisions between different value
    shapes (e.g., a list and a tuple with the same numeric contents).
    """
    if v is None or isinstance(v, (bool, int, str)):
        return v
    if isinstance(v, float):
        # NaN never compares equal under ``__eq__``, so the eq=>hash
        # invariant is vacuously satisfied. Map all NaNs to a single
        # sentinel so they at least hash deterministically.
        return ("__nan__",) if math.isnan(v) else v
    if isinstance(v, list):
        return ("list", tuple(_canon(x) for x in v))
    if hasattr(v, "to_numpy"):  # AMReX SmallMatrix (LinearMap, SpinMap)
        a = v.to_numpy()
        return ("ndarray", tuple(a.shape), tuple(a.ravel().tolist()))
    return ("repr", repr(v))


def _values_equal(a, b):
    if hasattr(a, "to_numpy") and hasattr(b, "to_numpy"):
        import numpy as np

        return np.array_equal(a.to_numpy(), b.to_numpy())
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_values_equal(x, y) for x, y in zip(a, b))
    return a == b


def _scalar_close(a, b, rtol, atol):
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) or math.isnan(b):
            return False
        return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)
    return a == b


def _values_close(a, b, rtol, atol):
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return False
        return all(_scalar_close(x, y, rtol, atol) for x, y in zip(a, b))
    if hasattr(a, "to_numpy") and hasattr(b, "to_numpy"):
        import numpy as np

        return np.allclose(a.to_numpy(), b.to_numpy(), rtol=rtol, atol=atol)
    return _scalar_close(a, b, rtol, atol)


def _dict_eq(d1, d2):
    if d1.keys() != d2.keys():
        return False
    return all(_values_equal(d1[k], d2[k]) for k in d1)


def _dict_close(d1, d2, rtol, atol):
    if d1.keys() != d2.keys():
        return False
    return all(_values_close(d1[k], d2[k], rtol, atol) for k in d1)


def _element_eq(self, other):
    """Value-based equality.

    Two elements are equal iff they are instances of the same class and
    their ``to_dict()`` outputs match key-for-key. Float values use
    ``==`` (so ``NaN != NaN``); list values are compared element-wise;
    AMReX SmallMatrix values use ``numpy.array_equal``.

    Returns ``NotImplemented`` when ``other`` is not the same element
    type, so Python's reflected-equality fallback applies (e.g., a
    foreign type's ``__eq__`` gets a chance, ultimately falling back to
    identity which yields ``False``).
    """
    if type(self) is not type(other):
        return NotImplemented
    return _dict_eq(self.to_dict(), other.to_dict())


def _element_hash(self):
    """Value-based hash, consistent with ``__eq__``.

    Two elements that compare equal under ``==`` produce the same hash.
    The hash reflects the element's *current* parameter values; mutating
    an element after using it as a ``set`` member or ``dict`` key
    invalidates the container, the same contract a hashable mutable
    object would have.
    """
    d = self.to_dict()
    return hash(
        (
            type(self).__name__,
            tuple(sorted((k, _canon(v)) for k, v in d.items())),
        )
    )


def _normalize_ignore(ignore_attributes):
    """Coerce ``ignore_attributes`` to a frozenset of key names.

    Accepts ``None``, a single string, or any iterable of strings.
    """
    if ignore_attributes is None:
        return frozenset()
    if isinstance(ignore_attributes, str):
        return frozenset({ignore_attributes})
    return frozenset(ignore_attributes)


def _element_isclose(
    self,
    other,
    *,
    rtol=_DEFAULT_RTOL,
    atol=_DEFAULT_ATOL,
    ignore_attributes=None,
):
    """Tolerant equality for lattice elements.

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
    ignore = _normalize_ignore(ignore_attributes)
    if "type" in ignore:
        if not hasattr(other, "to_dict"):
            return False
    elif type(self) is not type(other):
        return False
    d1 = {k: v for k, v in self.to_dict().items() if k not in ignore}
    d2 = {k: v for k, v in other.to_dict().items() if k not in ignore}
    return _dict_close(d1, d2, rtol, atol)


def isclose(a, b, **kwargs):
    """Free-function form of ``a.isclose(b)``.

    Provided so the call site reads symmetrically in ``a`` and ``b``,
    like :py:func:`math.isclose` / :py:func:`numpy.isclose`. Equivalent
    to the method form; forwards ``rtol``, ``atol``, and
    ``ignore_attributes``. Dispatches to ``b`` if ``a`` lacks
    ``isclose`` (e.g., a plain Python ``list`` paired with a
    ``KnownElementsList``).
    """
    if hasattr(a, "isclose"):
        return a.isclose(b, **kwargs)
    return b.isclose(a, **kwargs)


def register_elements_value_semantics(elements_module):
    """Attach ``__eq__``, ``__hash__``, and ``isclose`` to every class
    in ``elements_module`` that exposes a ``to_dict`` method, and the
    free-function ``isclose`` on the module itself.

    Parameters
    ----------
    elements_module : module
        Typically ``impactx.impactx_pybind.elements``.
    """
    for name in dir(elements_module):
        if name in _SKIP or name.startswith("_"):
            continue
        cls = getattr(elements_module, name)
        if not inspect.isclass(cls):
            continue
        if not hasattr(cls, "to_dict"):
            continue
        cls.__eq__ = _element_eq
        # Assigning ``__eq__`` from Python implicitly sets
        # ``__hash__ = None`` on the class, so we must reassign
        # ``__hash__`` afterwards to keep instances hashable.
        cls.__hash__ = _element_hash
        cls.isclose = _element_isclose
    elements_module.isclose = isclose
