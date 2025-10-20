"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell, Edoardo Zoni
License: BSD-3-Clause-LBNL
"""

import os
import re

from impactx import elements, kahan_babushka_sum


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
        from pals.BeamLine import BeamLine

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
    from pals.Drift import Drift
    from pals.Quadrupole import Quadrupole

    # Loop over the pals_beamline and create a new ImpactX KnownElementsList from it.
    #       Use self.extend(...) on the latter.
    ix_beamline = []
    for pals_element in pals_beamline.line:
        if isinstance(pals_element, Drift):
            ix_beamline.append(
                elements.Drift(
                    name=pals_element.name, ds=pals_element.length, nslice=nslice
                )
            )
        elif isinstance(pals_element, Quadrupole):
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


class FilteredElementsList:
    """A selection result class for ElementsList that maintains references to original elements.

    References to the original elements in a lattice are needed to allow modification of the original elements.
    """

    def __init__(self, original_list, indices):
        self._original_list = original_list
        self._indices = indices

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._original_list[self._indices[key]]
        elif isinstance(key, slice):
            # Return a new FilteredElementsList with sliced indices
            sliced_indices = self._indices[key]
            return FilteredElementsList(self._original_list, sliced_indices)
        else:
            raise TypeError(f"Invalid key type: {type(key)}")

    def __len__(self):
        return len(self._indices)

    def __iter__(self):
        for i in self._indices:
            yield self._original_list[i]

    def select(
        self,
        *,
        kind=None,
        name=None,
        s=None,
    ):
        """Apply filtering to this filtered list.

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
                     Examples: "quad1", r"quad\d+", ["quad1", "quad2"], [r"quad\d+", "bend1"]
        :type name: str or list[str] or tuple[str, ...] or None, optional

        :param s: Position range to filter by. Tuple/list with (lower, upper) bounds where None represents open-ended.
                  Elements are selected if ANY part overlaps with the range. Examples: (1.0, 5.0) for range 1.0 <= s <= 5.0, (1.0, None) for s >= 1.0, (None, 5.0) for s <= 5.0
        :type s: tuple[float | None, float | None] or list[float | None] or None, optional

        :return: FilteredElementsList containing references to original elements
        :rtype: FilteredElementsList

        :raises TypeError: If kind/name/s parameters have wrong types

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
                name=r"quad\d+"
            )  # Filter quads by regex pattern

            # Position-based filtering (always calculated from original lattice)
            early_elements = lattice.select(s=(None, 2.0))  # Elements with s <= 2.0
            drift_then_early = lattice.select(kind="Drift").select(
                s=(1.0, None)
            )  # Drift elements with s >= 1.0
        """
        # Apply filtering directly to the indices we already have
        if kind is not None or name is not None or s is not None:
            # Validate parameters
            _validate_select_parameters(kind, name, s)

            matching_indices = []

            for original_idx in self._indices:
                element = self._original_list[original_idx]
                if _check_element_match(
                    element, kind, name, s, original_idx, self._original_list
                ):
                    matching_indices.append(original_idx)

            return FilteredElementsList(self._original_list, matching_indices)

        # If no filtering criteria provided, return all current elements
        return FilteredElementsList(self._original_list, self._indices)

    def get_kinds(self) -> list[type]:
        """Get all unique element kinds in the filtered list.

        Returns:
            list[type]: List of unique element types (sorted by name).
        """
        return get_kinds(self)

    def count_by_kind(self, kind_pattern) -> int:
        """Count elements of a specific kind in the filtered list.

        Args:
            kind_pattern: The element kind to count. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            int: Number of elements of the specified kind.
        """
        return count_by_kind(self, kind_pattern)

    def has_kind(self, kind_pattern) -> bool:
        """Check if filtered list contains elements of a specific kind.

        Args:
            kind_pattern: The element kind to check for. Can be:
                - String name (e.g., "Drift", "Quad") - supports exact match
                - Regex pattern (e.g., r".*Quad") - supports pattern matching
                - Element type (e.g., elements.Drift) - supports exact type match

        Returns:
            bool: True if at least one element of the specified kind exists.
        """
        return has_kind(self, kind_pattern)

    def __repr__(self):
        return f"FilteredElementsList({len(self)} elements)"

    def __str__(self):
        return f"FilteredElementsList({len(self)} elements)"


def _is_regex_pattern(pattern: str) -> bool:
    """Check if a string looks like a regex pattern by testing if it contains regex metacharacters."""
    # Simple heuristic: if it contains regex metacharacters, treat as regex
    regex_chars = r".*+?^${}[]|()\\"
    return any(char in pattern for char in regex_chars)


def _matches_string(text: str, string_pattern: str) -> bool:
    """Check if text matches a string pattern (exact match or regex)."""
    if _is_regex_pattern(string_pattern):
        try:
            return bool(re.search(string_pattern, text))
        except re.error:
            # If regex compilation fails, fall back to exact match
            return text == string_pattern
    else:
        return text == string_pattern


def _validate_select_parameters(kind, name, s):
    """Validate parameters for select methods.

    Args:
        kind: Element type(s) to filter by
        name: Element name(s) to filter by
        s: Position range to filter by

    Raises:
        TypeError: If parameters have wrong types
    """
    if kind is not None:
        if isinstance(kind, (list, tuple)):
            for k in kind:
                if not isinstance(k, (str, type)):
                    raise TypeError(
                        "'kind' parameter must be a string/type or list of strings/types"
                    )
        elif not isinstance(kind, (str, type)):
            raise TypeError("'kind' parameter must be a string or type")

    if name is not None:
        if isinstance(name, (list, tuple)):
            for n in name:
                if not isinstance(n, str):
                    raise TypeError("'name' parameter must contain strings")
        elif not isinstance(name, str):
            raise TypeError(
                "'name' parameter must be a string or list/tuple of strings"
            )

    if s is not None:
        if isinstance(s, (list, tuple)):
            if len(s) != 2:
                raise TypeError(
                    "'s' parameter must have exactly 2 elements (lower, upper)"
                )
            for bound in s:
                if bound is not None and not isinstance(bound, (int, float)):
                    raise TypeError("'s' parameter bounds must be numbers or None")
        else:
            raise TypeError(
                "'s' parameter must be a tuple/list with 2 elements (lower, upper)"
            )


def _matches_kind_pattern(element, kind_pattern):
    """Check if an element matches a kind pattern.

    Args:
        element: The element to check
        kind_pattern: String or type to match against

    Returns:
        bool: True if element matches the pattern
    """
    if isinstance(kind_pattern, str):
        # String pattern (exact match or regex)
        return _matches_string(type(element).__name__, kind_pattern)
    elif isinstance(kind_pattern, type):
        # Element type (exact match)
        return type(element) is kind_pattern
    return False


def _matches_name_pattern(element, name_pattern):
    """Check if an element matches a name pattern.

    Args:
        element: The element to check
        name_pattern: String pattern to match against

    Returns:
        bool: True if element matches the pattern
    """
    return (
        hasattr(element, "has_name")
        and element.has_name
        and _matches_string(element.name, name_pattern)
    )


def _matches_s_position(element, s_range, element_s):
    """Check if an element's integrated position matches the s range criteria.

    An element matches if ANY part of it overlaps with the specified range.
    The element spans from element_s to element_s + element.ds.

    Args:
        element: The element to check
        s_range: Tuple/list of (lower, upper) bounds. None represents open-ended.
        element_s: The cumulative position of the element (calculated externally)

    Returns:
        bool: True if element's position overlaps with the range
    """
    if s_range is None:
        return True

    # Convert to tuple if it's a list
    if isinstance(s_range, list):
        s_range = tuple(s_range)

    if not isinstance(s_range, tuple) or len(s_range) != 2:
        raise TypeError(
            "'s' parameter must be a tuple/list with 2 elements (lower, upper)"
        )

    lower, upper = s_range

    # Element spans from element_s to element_s + element.ds
    element_start = element_s
    element_end = element_s + element.ds

    # Check if any part of the element overlaps with the range
    if lower is not None and element_end < lower:
        return False
    if upper is not None and element_start > upper:
        return False

    return True


def _check_element_match(element, kind, name, s, element_index, lattice):
    """Check if an element matches the given kind, name, and s criteria.

    Args:
        element: The element to check
        kind: Kind criteria (str, type, list, tuple, or None)
        name: Name criteria (str, list, tuple, or None)
        s: Position criteria (tuple/list with 2 elements, or None)
        element_index: Index of the element in the lattice
        lattice: The full lattice to calculate cumulative positions

    Returns:
        bool: True if element matches any criteria (OR logic)
    """
    match = False

    # Check for 'kind' parameter
    if kind is not None:
        if isinstance(kind, (str, type)):
            # Single kind
            if _matches_kind_pattern(element, kind):
                match = True
        elif isinstance(kind, (list, tuple)):
            # Multiple kinds (OR logic)
            for k in kind:
                if _matches_kind_pattern(element, k):
                    match = True
                    break

    # Check for 'name' parameter (only if kind didn't match - OR logic)
    if name is not None and not match:
        if isinstance(name, str):
            # Single name
            if _matches_name_pattern(element, name):
                match = True
        elif isinstance(name, (list, tuple)):
            # Multiple names (OR logic)
            for name_item in name:
                if _matches_name_pattern(element, name_item):
                    match = True
                    break

    # Check for 's' parameter (only if neither kind nor name matched - OR logic)
    if s is not None and not match:
        # Calculate cumulative position up to this element using accurate summation
        ds_values = [lattice[i].ds for i in range(element_index)]
        cumulative_s = kahan_babushka_sum(ds_values)

        if _matches_s_position(element, s, cumulative_s):
            match = True

    return match


def select(
    self,
    *,
    kind=None,
    name=None,
    s=None,
) -> FilteredElementsList:
    """Filter elements by type, name, and position with OR-based logic.

    This method supports filtering elements by their type, name, and/or integrated position using keyword arguments.
    Returns references to original elements, allowing modification and chaining.

    **Filtering Logic:**

    - **Within a single filter**: OR logic (e.g., ``kind=["Drift", "Quad"]`` matches Drift OR Quad)
    - **Between different filters**: OR logic (e.g., ``kind="Quad", name="quad1", s=(1.0, 5.0)`` matches Quad OR named "quad1" OR in position range)
    - **Chaining filters**: AND logic (e.g., ``lattice.select(kind="Drift").select(name="drift1")`` matches Drift AND named "drift1")

    :param kind: Element type(s) to filter by. Can be a single string/type or a list/tuple
                 of strings/types for OR-based filtering. String values support exact matches
                 and regex patterns. Examples: "Drift", r".*Quad", elements.Drift, ["Drift", r".*Quad"], [elements.Drift, elements.Quad]
    :type kind: str or type or list[str | type] or tuple[str | type, ...] or None, optional

    :param name: Element name(s) to filter by. Can be a single string, regex pattern string, or
                 a list/tuple of strings and/or regex pattern strings for OR-based filtering.
                 Examples: "quad1", r"quad\d+", ["quad1", "quad2"], [r"quad\d+", "bend1"]
    :type name: str or list[str] or tuple[str, ...] or None, optional

    :param s: Position range to filter by. Tuple/list with (lower, upper) bounds where None represents open-ended.
              Elements are selected if ANY part overlaps with the range. Examples: (1.0, 5.0) for range 1.0 <= s <= 5.0, (1.0, None) for s >= 1.0, (None, 5.0) for s <= 5.0
    :type s: tuple[float | None, float | None] or list[float | None] or None, optional

    :return: FilteredElementsList containing references to original elements
    :rtype: FilteredElementsList

    :raises TypeError: If kind/name/s parameters have wrong types

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

         lattice.select(name=r"quad\d+")  # Get elements matching pattern
         lattice.select(name=[r"quad\d+", "bend1"])  # Mix regex and strings

     Position-based filtering:

     .. code-block:: python

         lattice.select(s=(1.0, 5.0))  # Elements that overlap with range 1.0 <= s <= 5.0
         lattice.select(s=(1.0, None))  # Elements that overlap with s >= 1.0
         lattice.select(s=(None, 5.0))  # Elements that overlap with s <= 5.0
         lattice.select(kind="Drift", s=(0.0, 2.0))  # Drift elements OR overlapping range 0.0 <= s <= 2.0

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

    # Handle keyword arguments for filtering
    if kind is not None or name is not None or s is not None:
        # Validate parameters
        _validate_select_parameters(kind, name, s)

        matching_indices = []

        for i, element in enumerate(self):
            if _check_element_match(element, kind, name, s, i, self):
                matching_indices.append(i)

        return FilteredElementsList(self, matching_indices)

    # If no filtering criteria provided, return all elements
    all_indices = list(range(len(self)))
    return FilteredElementsList(self, all_indices)


def get_kinds(self) -> list[type]:
    """Get all unique element kinds in the list.

    Returns:
        list[type]: List of unique element types (sorted by name).
    """
    kinds = set()
    for element in self:
        kinds.add(type(element))
    return sorted(list(kinds), key=lambda t: t.__name__)


def count_by_kind(self, kind_pattern) -> int:
    """Count elements of a specific kind.

    Args:
        kind_pattern: The element kind to count. Can be:
            - String name (e.g., "Drift", "Quad") - supports exact match
            - Regex pattern (e.g., r".*Quad") - supports pattern matching
            - Element type (e.g., elements.Drift) - supports exact type match

    Returns:
        int: Number of elements of the specified kind.
    """
    count = 0
    for element in self:
        if _matches_kind_pattern(element, kind_pattern):
            count += 1
    return count


def has_kind(self, kind_pattern) -> bool:
    """Check if list contains elements of a specific kind.

    Args:
        kind_pattern: The element kind to check for. Can be:
            - String name (e.g., "Drift", "Quad") - supports exact match
            - Regex pattern (e.g., r".*Quad") - supports pattern matching
            - Element type (e.g., elements.Drift) - supports exact type match

    Returns:
        bool: True if at least one element of the specified kind exists.
    """
    for element in self:
        if _matches_kind_pattern(element, kind_pattern):
            return True
    return False


def register_KnownElementsList_extension(kel):
    """KnownElementsList helper methods"""
    from ..plot.Survey import plot_survey

    # register member functions for KnownElementsList
    kel.from_pals = from_pals
    kel.load_file = load_file
    kel.plot_survey = plot_survey

    # Enhanced element selection methods
    kel.select = select
    kel.get_kinds = get_kinds
    kel.count_by_kind = count_by_kind
    kel.has_kind = has_kind
