"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import inspect
from typing import Dict, List, Type


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

        return docstrings
