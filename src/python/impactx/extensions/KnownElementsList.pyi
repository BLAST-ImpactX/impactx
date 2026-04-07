"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell, Edoardo Zoni
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

import math as math
import os as os
import re as re
import weakref as weakref

import impactx.impactx_pybind.elements
from impactx.impactx_pybind import elements
from impactx.impactx_pybind.elements import FilteredElementsList

__all__: list[str] = [
    "FILTERED_ELEMENTS_LIST_INVALID_MSG",
    "FilteredElementsList",
    "count_by_kind",
    "elements",
    "from_dicts",
    "from_pals",
    "get_kinds",
    "has_kind",
    "load_file",
    "math",
    "os",
    "re",
    "register_KnownElementsList_extension",
    "select",
    "to_dicts",
    "to_py",
    "weakref",
]

def _check_element_match(element, kind, name):
    """
    Check if an element matches the given kind and name criteria.

    Args:
        element: The element to check
        kind: Kind criteria (str, type, list, tuple, or None)
        name: Name criteria (str, list, tuple, or None)

    Returns:
        bool: True if element matches any criteria (OR logic)
    """

def _clone_element(template):
    """
    Deep-clone a lattice element via ``to_dict`` (pybind elements are not copy.copy-able).
    """

def _commit_lattice_rebuild(original, new_elements) -> None:
    """
    Replace lattice contents with ``new_elements`` and invalidate all FilteredElementsList views.
    """

def _drift_class_for_replace_with_drifts(model: str, old_el) -> type:
    """
    Map ``model`` and ``old_el`` to the Drift / ChrDrift / ExactDrift class to insert.

    For ``model=="match"``, the class follows ``_model_key_from_element_typename``; otherwise
    ``model`` must already be validated against ``_DRIFT_MODEL_CLASSES``.
    """

def _element_from_dict(d: dict):
    """
    Create an element from its to_dict() representation.

    Args:
        d: Dictionary from element.to_dict(), must include 'type' key

    Returns:
        Element instance of the appropriate type

    Raises:
        KeyError: If 'type' key is missing
        AttributeError: If element type is not found in elements module
    """

def _filter_kwargs(d: dict) -> dict:
    """
    Filter a to_dict() result into valid constructor kwargs.

    Removes 'type' (not a constructor argument) and 'ds' when zero
    (thin elements don't accept ds).

    Args:
        d: Dictionary from element.to_dict()

    Returns:
        dict: Filtered dictionary suitable for element constructor
    """

def _format_value(v):
    """
    Format a value for Python code generation.

    Converts AMReX SmallMatrix objects to lists.
    Returns a tuple of (formatted_value, matrix_dims) where matrix_dims
    is (rows, cols) from the SmallMatrix type, or None otherwise.
    """

def _invalidate_all_registered_views(lattice) -> None:
    """
    Mark every registered FilteredElementsList for this lattice as invalid.
    """

def _is_regex_pattern(pattern: str) -> bool:
    """
    Check if a string looks like a regex pattern by testing if it contains regex metacharacters.
    """

def _make_drift_from_old(
    cls, old_el, *, keep_name, keep_ds, keep_alignment, keep_aperture
):
    """
    Build a drift (``cls`` is Drift / ChrDrift / ExactDrift) from thick-element fields on ``old_el``.

    When ``keep_ds`` is False, ``ds`` is set to 0. When ``keep_name`` is False, ``name`` is None.
    If ``keep_ds`` is True but ``old_el`` has no ``ds`` attribute (thin element), ``ds`` defaults to 0.
    When ``keep_alignment`` is True, copy dx/dy/rotation from the old element.
    When ``keep_aperture`` is True, copy aperture_x/aperture_y from the old element.
    """

def _matches_kind_pattern(element, kind_pattern):
    """
    Check if an element matches a kind pattern.

    Args:
        element: The element to check
        kind_pattern: String or type to match against

    Returns:
        bool: True if element matches the pattern
    """

def _matches_name_pattern(element, name_pattern):
    """
    Check if an element matches a name pattern.

    Args:
        element: The element to check
        name_pattern: String pattern to match against

    Returns:
        bool: True if element matches the pattern
    """

def _matches_string(text: str, string_pattern: str) -> bool:
    """
    Check if text matches a string pattern (exact match or regex).
    """

def _model_key_from_element_typename(type_name: str) -> str:
    """
    Return the drift-model key for an element class name (linear / paraxial / exact).
    """

def _rad2deg(radians: float) -> float:
    """
    Convert radians to degrees.
    """

def _registry_for(lattice):
    """
    Return the WeakSet of FilteredElementsList instances for this lattice.
    """

def _validate_select_parameters(kind, name):
    """
    Validate parameters for select methods.

    Args:
        kind: Element type(s) to filter by
        name: Element name(s) to filter by

    Raises:
        TypeError: If parameters have wrong types
    """

def count_by_kind(self, kind_pattern) -> int:
    """
    Count elements of a specific kind.

    Args:
        kind_pattern: The element kind to count. Can be:
            - String name (e.g., "Drift", "Quad") - supports exact match
            - Regex pattern (e.g., r".*Quad") - supports pattern matching
            - Element type (e.g., elements.Drift) - supports exact type match

    Returns:
        int: Number of elements of the specified kind.
    """

def from_dicts(self, dicts: list[dict]):
    """
    Load and append elements from a list of dictionaries.

    Each dictionary should be in the format produced by element.to_dict(),
    containing at minimum a 'type' key identifying the element class.

    Args:
        dicts: List of element dictionaries

    Example:
        .. code-block:: python

            import json
            from impactx import elements

            # Load from JSON
            with open("lattice.impactx.json") as f:
                data = json.load(f)

            lattice = elements.KnownElementsList()
            lattice.from_dicts(data)

    Note:
        Elements with matrix parameters (LinearMap, SpinMap) require
        the matrices to be AMReX SmallMatrix objects. Use
        :func:`impactx.extensions.matrix_hook` as a JSON object_hook
        when loading such elements.
    """

def from_pals(self, pals_beamline, nslice=1):
    """
    Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) Python BeamLine.

    https://github.com/campa-consortium/pals-python
    """

def get_kinds(self) -> list[type]:
    """
    Get all unique element kinds in the list.

    Returns:
        list[type]: List of unique element types (sorted by name).
    """

def has_kind(self, kind_pattern) -> bool:
    """
    Check if list contains elements of a specific kind.

    Args:
        kind_pattern: The element kind to check for. Can be:
            - String name (e.g., "Drift", "Quad") - supports exact match
            - Regex pattern (e.g., r".*Quad") - supports pattern matching
            - Element type (e.g., elements.Drift) - supports exact type match

    Returns:
        bool: True if at least one element of the specified kind exists.
    """

def load_file(self, filename, nslice=1):
    """
    Load and append a lattice file from MAD-X (.madx) or PALS (e.g., .pals.yaml) formats.
    """

def register_KnownElementsList_extension(kel):
    """
    KnownElementsList helper methods
    """

def select(
    self, *, kind=None, name=None
) -> impactx.impactx_pybind.elements.FilteredElementsList:
    """
    Filter elements by type and name with OR-based logic.

    This method supports filtering elements by their type and/or name using keyword arguments.
    Returns references to original elements, allowing modification and chaining.

    **Filtering Logic:**

    - **Within a single filter**: OR logic (e.g., ``kind=["Drift", "Quad"]`` matches Drift OR Quad)
    - **Between different filters**: OR logic (e.g., ``kind="Quad", name="quad1"`` matches Quad OR named "quad1")
    - **Chaining filters**: AND logic (e.g., ``lattice.select(kind="Drift").select(name="drift1")`` matches Drift AND named "drift1")

    :param kind: Element type(s) to filter by. Can be a single string/type or a list/tuple
                 of strings/types for OR-based filtering. String values support exact matches
                 and regex patterns. Examples: "Drift", r".*Quad", elements.Drift, ["Drift", r".*Quad"], [elements.Drift, elements.Quad]
    :type kind: str or type or list[str | type] or tuple[str | type, ...] or None, optional

    :param name: Element name(s) to filter by. Can be a single string, regex pattern string, or
                 a list/tuple of strings and/or regex pattern strings for OR-based filtering.
                 Examples: "quad1", r"quad\\d+", ["quad1", "quad2"], [r"quad\\d+", "bend1"]
    :type name: str or list[str] or tuple[str, ...] or None, optional

    :return: FilteredElementsList containing references to original elements
    :rtype: FilteredElementsList

    :raises TypeError: If kind/name parameters have wrong types

    **Examples:**

    Single value filtering:

    .. code-block:: python

        lattice.select(kind="Drift")  # Get all drift elements (string)
        lattice.select(kind=elements.Drift)  # Get all drift elements (type)
        lattice.select(
            kind=r".*Quad"
        )  # Get all elements matching regex pattern (Quad, ExactQuad, ChrQuad)
        lattice.select(name="quad1")  # Get elements named "quad1"
        lattice.select(
            kind="Quad", name="quad1"
        )  # Get quad elements OR elements named "quad1"

    OR-based filtering with lists (within single filter):

    .. code-block:: python

        lattice.select(kind=["Drift", "Quad"])  # Get drift OR quad elements (strings)
        lattice.select(kind=[elements.Drift, elements.Quad])  # Get drift OR quad elements (types)
        lattice.select(kind=["Drift", elements.Quad])  # Mix strings and types
        lattice.select(kind=[r".*Quad", r".*Bend.*"])  # Mix regex patterns
        lattice.select(name=["quad1", "quad2"])  # Get elements named "quad1" OR "quad2"

     Regex pattern filtering:

     .. code-block:: python

         lattice.select(name=r"quad\\d+")  # Get elements matching pattern
         lattice.select(name=[r"quad\\d+", "bend1"])  # Mix regex and strings

    Chaining filters (AND logic between chained calls):

    .. code-block:: python

        lattice.select(kind="Drift").select(
            name="drift1"
        )  # Drift elements AND named "drift1"
        lattice.select(kind="Quad")[0]  # First quad element
        lattice.select(name="quad1").select(
            kind="Quad"
        )  # Elements named "quad1" AND of type "Quad"

    Reference preservation and modification:

    .. code-block:: python

        drift_elements = lattice.select(kind="Drift")
        drift_elements[0].ds = 5.0  # Modifies the original element in lattice
        assert lattice[0].ds == 5.0  # Original element is modified

    Modification of elements (reference preservation):

    .. code-block:: python

        drift = lattice.select(kind="Drift")[0]  # Get first drift element
        drift.ds = 2.0  # Modify original element
        quad_elements = lattice.select(kind="Quad")  # Get all quad elements
        quad_elements[0].k = 1.5  # Modify first quad's strength
        # All modifications affect the original lattice elements
    """

def to_dicts(self) -> list[dict]:
    """
    Serialize the lattice to a list of dictionaries.

    Each element is converted to a dictionary using its to_dict() method.
    The resulting list can be serialized to JSON, YAML, or other formats.

    Returns:
        list[dict]: List of element dictionaries

    Example:
        .. code-block:: python

            import json
            from impactx import elements

            lattice = elements.KnownElementsList(
                [
                    elements.Drift(ds=1.0, name="d1"),
                    elements.Quad(ds=0.5, k=2.0, name="q1"),
                ]
            )

            # Serialize to JSON
            data = lattice.to_dicts()
            with open("lattice.impactx.json", "w") as f:
                json.dump(data, f, indent=2)

    Note:
        Elements with matrix parameters (LinearMap, SpinMap) contain
        AMReX SmallMatrix objects that require custom JSON encoding.
        Use :func:`impactx.extensions.ImpactXEncoder` for JSON serialization
        of such elements.
    """

def to_py(self) -> str:
    """
    Generate Python code that recreates this lattice.

    Returns a string containing a complete Python script with imports
    and a ``get_lattice()`` function that returns a KnownElementsList
    with all elements.

    Returns:
        str: Python source code

    Example:
        .. code-block:: python

            from impactx import elements

            lattice = elements.KnownElementsList(
                [
                    elements.Drift(ds=1.0, name="d1"),
                    elements.Quad(ds=0.5, k=2.0, name="q1"),
                ]
            )

            # Generate Python code
            code = lattice.to_py()
            print(code)

            # Save to file
            with open("my_lattice.py", "w") as f:
                f.write(code)

            # Later, use the generated file:
            # from my_lattice import get_lattice
            # lattice = get_lattice()
    """

FILTERED_ELEMENTS_LIST_INVALID_MSG: str = "This lattice selection is no longer valid because the lattice was modified; call select() again on the lattice."
_DRIFT_MODEL_CLASSES: dict = {
    "linear": impactx.impactx_pybind.elements.Drift,
    "paraxial": impactx.impactx_pybind.elements.ChrDrift,
    "exact": impactx.impactx_pybind.elements.ExactDrift,
}
_filtered_views_by_lattice: weakref.WeakKeyDictionary  # value = <WeakKeyDictionary>
