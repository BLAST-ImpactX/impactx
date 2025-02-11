import inspect


class InputDefaultsHelper:
    """
    Methods in this class are used to dynamically parse
    core ImpactX data (default values, docstrings, etc.)
    """

    @staticmethod
    def get_docstrings(classes, default_list) -> dict:
        """
        Retrieves docstrings for each method and property
        in the provided clases.

        :param classes: The class names to parse docstrings with.
        :param defaults_list: The dictionary of defaults value.
        """

        docstrings = {}

        for each_class in classes:
            for name, attribute in inspect.getmembers(each_class):
                if name not in default_list:
                    continue

                is_method = inspect.isfunction(attribute)
                is_property = inspect.isdatadescriptor(attribute)

                if is_method or is_property:
                    docstring = inspect.getdoc(attribute) or ""
                    docstrings[name] = docstring

        return docstrings
