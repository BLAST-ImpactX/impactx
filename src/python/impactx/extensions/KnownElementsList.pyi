"""

This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell, Edoardo Zoni
License: BSD-3-Clause-LBNL
"""

from __future__ import annotations

import os as os
import re as re

from impactx.impactx_pybind import elements

__all__: list[str] = [
    "FilteredElementsList",
    "count_by_kind",
    "elements",
    "from_pals",
    "get_kinds",
    "has_kind",
    "load_file",
    "os",
    "re",
    "register_KnownElementsList_extension",
    "select",
]

class FilteredElementsList:
    """
    A selection result class for ElementsList that maintains references to original elements.

    References to the original elements in a lattice are needed to allow modification of the original elements.
    """
    def __getitem__(self, key): ...
    def __init__(self, original_list, indices): ...
    def __iter__(self): ...
    def __len__(self): ...
    def __repr__(self): ...
    def __str__(self): ...
    def count_by_kind(self, kind_pattern) -> int:
        """
        Count elements of a specific kind in the filtered list.

        Args:
            kind_pattern: The element kind to count. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            int: Number of elements of the specified kind.
        """
    def get_kinds(self) -> list[type]:
        """
        Get all unique element kinds in the filtered list.

        Returns:
            list[type]: List of unique element types (sorted by name).
        """
    def has_kind(self, kind_pattern) -> bool:
        """
        Check if filtered list contains elements of a specific kind.

        Args:
            kind_pattern: The element kind to check for. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            bool: True if at least one element of the specified kind exists.
        """
    def select(self, *, kind=None, name=None):
        """
        Apply filtering to this filtered list.

        This method applies additional filtering to an already filtered list,
        maintaining references to the original elements and enabling chaining.

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

        Additional filtering on already filtered results:

        .. code-block:: python

            drift_elements = lattice.select(
                kind="Drift"
            )  # or lattice.select(kind=elements.Drift)
            first_drift = drift_elements.select(
                name="drift1"
            )  # Further filter drifts by name
            quad_elements = lattice.select(
                kind="Quad"
            )  # or lattice.select(kind=elements.Quad)
            strong_quads = quad_elements.select(
                name=r"quad\\d+"
            )  # Filter quads by regex pattern
        """

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

def _is_regex_pattern(pattern: str) -> bool:
    """
    Check if a string looks like a regex pattern by testing if it contains regex metacharacters.
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

def select(self, *, kind=None, name=None) -> FilteredElementsList:
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
