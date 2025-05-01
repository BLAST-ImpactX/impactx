"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import inspect
import re
from typing import Callable, Dict, List, Type

from impactx.distribution_input_helpers import twiss


class InputDefaultsHelper:
    """
    Methods in this class are used to dynamically parse
    core ImpactX data (default values, docstrings, etc.)
    """

    @staticmethod
    def get_docstrings(
        class_names: List[Type], default_list: Dict[str, any]
    ) -> Dict[str, str]:
        """
        Retrieves docstrings for each method and property
        in the provided clases.

        :param classes: The class names to parse docstrings with.
        :param defaults_list: The dictionary of defaults value.
        """

        docstrings = {}

        for each_class in class_names:
            for name, attribute in inspect.getmembers(each_class):
                if name not in default_list:
                    continue

                is_method = inspect.isfunction(attribute)
                is_property = inspect.isdatadescriptor(attribute)

                if is_method or is_property:
                    docstring = inspect.getdoc(attribute) or ""
                    docstrings[name] = docstring

        distribution_tooltips = InputDefaultsHelper.get_tooltips_from_param(twiss)
        docstrings.update(distribution_tooltips)

        return docstrings

    @staticmethod
    def get_tooltips_from_param(func: Callable) -> Dict[str, str]:
        """
        Extract all ':param name: description' entries from a function's docstring.

        :param func: The function whose docstring you want to parse.
        :return: A dict mapping each parameter name to its description.
        """
        tooltips = {}
        doc = inspect.getdoc(func) or ""
        pattern = re.compile(r"^\s*:param\s+(\w+)\s*:\s*(.+)$", re.MULTILINE)
        pattern_matches = list(pattern.finditer(doc))

        if not pattern_matches:
            raise ValueError(f"Found no docstrings to parse in function {func.__name__}")

        for match in pattern_matches:
            param_name = match.group(1)
            description = match.group(2)
            tooltips[param_name] = description

        return tooltips
