import re
from typing import Any, Dict, List


class DashboardLatticeParser:
    """
    Helper class to parse the lattice configuration from Impactx .py-compatible simulation files.
    """

    def __init__(self, content: str):
        self._flatten_cache: Dict[str, List[str]] = {}
        self._variable_assignments_cache: Dict[str, Dict[str, str]] = {}
        self._content = content
        self._content_hash: str = str(hash(content))

    def parse(self) -> Dict[str, Any]:
        """
        Extracts the lattice configuration from provided ImpactX simulation file.

        Example return:
        {
            "lattice_elements": [
                {"name": "Drift", "parameters": {"ds": "0.5"}},
                {"name": "Quad", "parameters": {"k": "quad_strength", "ds": "0.3"}}
            ]
        }

        :return: Parsed dictionary containing 'lattice_elements'.
        """
        lattice_order = self.collect_lattice_operations(debug=False)

        expanded_elements = []
        for operation in lattice_order:
            operation_type = operation["type"]
            operation_arg = operation["argument"]

            match operation_type:
                case "extend":
                    expanded_elements.extend(self._flatten(operation_arg))
                case "append":
                    expanded_elements.append(operation_arg)
                case "reverse":
                    self._apply_reverse(operation_arg)
                case _:
                    print(f"Warning: Unsupported operation type: {operation_type}")

        clean_lattice_list = self.replace_variables(expanded_elements)
        clean_lattice_list_str = "\n".join(clean_lattice_list)
        result = self.parse_cleaned_lattice(clean_lattice_list_str)

        return result

    def parse_cleaned_lattice(self, content: str) -> Dict[str, List[Dict[str, Any]]]:
        """
        Parses the lattice elements from the ImpactX simulation file content.

        Extracts element names and their parameters from constructor calls in the format:
        elements.ElementName(param1=value1, param2=value2, ...)

        EX:
            elements.Drift(ds=1.0)

            Results in:
            {
                "lattice_elements": [
                    {
                        "name": "Drift",
                        "parameters": {"ds": "1.0"}
                    }
                ]
            }

        :param content: The content of the ImpactX simulation file.
        :return: A dictionary containing the parsed lattice elements.
        """

        dictionary = {"lattice_elements": []}

        element_pattern = r"elements\.(\w+)\((.*?)\)"  # EX: elements.Drift(...)
        lattice_elements = re.findall(element_pattern, content, re.DOTALL)

        for element_name, parameter in lattice_elements:
            element = {"name": element_name, "parameters": {}}

            # CHANGE: Updated parameter pattern to handle multiline and whitespace around =
            # OLD: r"(\w+)=([^,\)]+)"
            # NEW: r"(\w+)\s*=\s*([^,\)\n]+)" with re.MULTILINE flag
            parameter_pattern = r"(\w+)\s*=\s*([^,\)\n]+)"  # EX: ds=1.0, k1=0.5
            all_parameters = re.findall(parameter_pattern, parameter)

            for parameter_name, value in all_parameters:
                element["parameters"][parameter_name] = value.strip("'\"")

            dictionary["lattice_elements"].append(element)

        return dictionary

    def collect_lattice_operations(self, debug: bool = False) -> List[Dict[str, str]]:
        """
        Collects lattice operations (sim.lattice.append(), sim.lattice.extend(), and variable.reverse() calls)
        in the order they appear in the content.

        EX:
            sim.lattice.append(monitor) ; sim.lattice.extend([drift1, quad1]) ; lattice_half.reverse()

            Results the following (in order):
            [
                {"type": "append", "argument": "monitor"},
                {"type": "extend", "argument": "[drift1, quad1]"}
                {"type": "reverse", "argument": "lattice_half"}
            ]

        :param debug: Whether to print the collected operations.
        :return: List of dictionaries with 'type' and 'argument' keys.
        """
        operations = []

        # Captures the operation type (append or extend) and its argument
        # (everything inside the parentheses, up to the last closing parenthesis on the line).
        lattice_call_pattern = r"sim\.lattice\.(append|extend)\((.*)\)"

        # Store sim.lattice.append and sim.lattice.extend calls
        for match in re.finditer(lattice_call_pattern, self._content):
            operation, arg = match.groups()
            operations.append(
                (match.start(), {"type": operation, "argument": arg.strip()})
            )

        # Store .reverse() calls
        reverse_pattern = r"(\w+)\.reverse\(\)"
        for match in re.finditer(reverse_pattern, self._content):
            operations.append(
                (match.start(), {"type": "reverse", "argument": match.group(1)})
            )

        # important: sort operations by their position in the content
        # since the for loops can be executed in any order
        operations = [type for _, type in sorted(operations, key=lambda x: x[0])]

        if debug:
            self.print_lattice_operations(operations)

        return operations

    def _get_variable_assignments(self) -> Dict[str, str]:
        """
        Helper function to extract all variable list assignments from content.
        Caches the result to avoid re-parsing the same content.

        EX:
            content = '''
                drift1 = elements.Drift(ds=1.0)
                cell = [drift1, quad1]
                line = [cell, monitor]
            '''
            _get_variable_assignments()

            Results in:
                {
                    "cell": "drift1, quad1",
                    "line": "cell, monitor"
                }

        """
        if self._content_hash in self._variable_assignments_cache:
            return self._variable_assignments_cache[self._content_hash]

        variable_assignments = {}
        var_assignment_pattern = r"(\w+)\s*=\s*\[(.*?)\]"

        for match in re.finditer(var_assignment_pattern, self._content, re.DOTALL):
            var_name = match.group(1)
            list_content = match.group(2)
            variable_assignments[var_name] = list_content

        self._variable_assignments_cache[self._content_hash] = variable_assignments
        return variable_assignments

    def _flatten(self, variable_name: str, debug: bool = False) -> List[str]:
        """
        Recursively expands a varible name to replace it with its set of elements.
        Utilizes caching to avoid redundant parsing.

        EX:
            content = '''
                cell = [drift1, quad1]
                line = [cell, cell]
                sim.lattice.extend(line)
            '''
            _flatten(content, "line")

            Results in:
                ["drift1", "quad1", "drift1", "quad1"]

        :param variable_name: Name of the specific variable to expand (e.g. "line")
        :param debug: Whether to print the expanded list.
        :return: List of individual element names with all nesting resolved.
        """
        # check cache first
        cache_key = f"{self._content_hash}:{variable_name}"
        if cache_key in self._flatten_cache:
            return self._flatten_cache[cache_key]

        # Get variable assignments (cached)
        variable_assignments = self._get_variable_assignments()

        # Check if the input is an inline list like "[monitor, elements.Drift(...)]"
        if variable_name.startswith("[") and variable_name.endswith("]"):
            list_contents = variable_name[1:-1]
            # split on commas that are NOT inside parentheses
            list_to_flatten = [
                element.strip()
                for element in re.split(r",\s*(?![^()]*\))", list_contents)
                if element.strip()
            ]
        else:
            if variable_name not in variable_assignments:
                self._flatten_cache[cache_key] = [variable_name]
                return [variable_name]  # It's not a list, it's a single element

            list_content = variable_assignments[variable_name]
            list_to_flatten = [
                item.strip() for item in list_content.split(",") if item.strip()
            ]

        expanded = []
        for item in list_to_flatten:
            # recursively expand each item
            sub_items = self._flatten(item, debug)
            expanded.extend(sub_items)

        # Cache the result
        self._flatten_cache[cache_key] = expanded

        return expanded

    def replace_variables(self, raw_lattice: List[str]) -> List[str]:
        """
        This function is called to simplify the lattice list by replacing variable names with their corresponding constructor calls.

        EX:
            (input)
                drift1 = elements.Drift(ds=1.0)
                quad1 = elements.Quad(k=0.5)
                raw_lattice = ["drift1", "quad1"]
            (output)
                raw_lattice = ["elements.Drift(ds=1.0)", "elements.Quad(k=0.5)"]

        :param raw_lattice: List of lattice element variable names or constructor calls, e.g. ["drift1", "quad1"].
        :return: List with variable names replaced by their corresponding constructor calls.
        """
        element_mapping = {}

        ellement_assignment_pattern = r"^\s*(\w+)\s*=\s*(elements\.\w+\(.*?\))"
        all_element_assignments = re.findall(
            ellement_assignment_pattern, self._content, re.MULTILINE | re.DOTALL
        )

        for var_name, element in all_element_assignments:
            element_mapping[var_name] = element

        if not element_mapping:
            return raw_lattice

        # Replace each item in the list
        # later can be optimized by not iterating over the whole raw_lattice list
        return [element_mapping.get(item, item) for item in raw_lattice]

    def extract_lattice_inputs(self, parsed_lattice: Dict[str, Any]) -> List[str]:
        """
        Extracts all parameter values from parsed lattice data.
        """
        used_variables = set()

        for element in parsed_lattice.get("lattice_elements", []):
            for value in element["parameters"].values():
                used_variables.add(str(value).strip())

        return list(used_variables)

    def _apply_reverse(self, var_name: str) -> None:
        """
        Reverse the order of elements inside a list variable (in place).

        Example:
            1.) line = [a, b, c]
            2.) line.reverse()
            3.) result is -> [c, b, a]

        :param var_name: Name of the variable to reverse.
        """

        variable_assignments = self._get_variable_assignments()
        if var_name not in variable_assignments:
            return

        # Take the string inside the brackets, e.g. "a, b, c".
        list_content = variable_assignments[var_name]

        # Split that string into individual items.
        # We split on commas, but the regex makes sure we do not split commas
        # that are inside parentheses (e.g., elements.Quad(k=0.5, ds=1.0)).
        raw_items = re.split(r",\s*(?![^()]*\))", list_content)

        items = [item.strip() for item in raw_items if item.strip()]
        items.reverse()
        variable_assignments[var_name] = ", ".join(items)

        # Clear the flatten cache so future lookups rebuild the list using the new reversed order.
        self._flatten_cache.clear()

    # -----------------------------------------------------------------------------
    # Debug methods
    # -----------------------------------------------------------------------------

    def print_lattice_operations(self, operations: List[Dict[str, str]]) -> None:
        """
        Prints all lattice operations in the order they appear.
        """

        print("Full lattice operation sequence (in order):")
        for operation in operations:
            print(f"  {operation['type']}({operation['argument']})")
