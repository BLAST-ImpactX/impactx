import re
from typing import Dict, List, Set, Any

class DashboardLatticeConfigParser:
    """
    Helper class to parse the lattice configuration from Impactx .py-compatible simulation files.
    """

    def parse_lattice(self, content: str) -> Dict[str, Any]:
        """
        Parses lattice elements and also extracts used lattice variables.
        """
        lattice_order = self.collect_lattice_operations(content, debug=True)

        expanded_elements = []
        for operation in lattice_order:
            operation_type = operation["type"]
            operation_arg = operation["argument"]

            match operation_type:
                case "extend":
                    expanded_elements.extend(self._flatten(content, operation_arg))
                case "append":
                    expanded_elements.append(operation_arg)
                case "reverse":
                    # Find the variable definition and reverse its flattened list
                    variable_elements = self._flatten(content, operation_arg)
                    expanded_elements.extend(reversed(variable_elements))
                case _:
                    print(f"Warning: Unsupported operation type: {operation_type}")

        clean_lattice_list = self.replace_variables(content, expanded_elements)
        clean_lattice_list_str = '\n'.join(clean_lattice_list)
        result = self.parse_cleaned_lattice(clean_lattice_list_str)
        
        result["used_lattice_variables"] = self.extract_used_variables(result)

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
            element = {
                "name": element_name,
                "parameters": {}
            }

            # CHANGE: Updated parameter pattern to handle multiline and whitespace around =
            # OLD: r"(\w+)=([^,\)]+)"
            # NEW: r"(\w+)\s*=\s*([^,\)\n]+)" with re.MULTILINE flag
            parameter_pattern = r"(\w+)\s*=\s*([^,\)\n]+)"  # EX: ds=1.0, k1=0.5
            all_parameters = re.findall(parameter_pattern, parameter)

            for parameter_name, value in all_parameters:
                element["parameters"][parameter_name] = value.strip("'\"")

            dictionary["lattice_elements"].append(element)

        return dictionary


    def collect_lattice_operations(self, content: str, debug: bool = False) -> List[Dict[str, str]]:
        """
        Returns a list of operations (in order) that define how the lattice is built.
        Handles sim.lattice.append, sim.lattice.extend, and .reverse() calls.

        EX:
            sim.lattice.append(monitor)
            sim.lattice.extend([drift1, quad1])
            lattice_half.reverse()

            Results the following (in order):
            [
                {"type": "append", "argument": "monitor"},
                {"type": "extend", "argument": "[drift1, quad1]"}
                {"type": "reverse", "argument": "lattice_half"}
            ]

        :param content: Full text content of the ImpactX simulation file.
        :param debug: Whether to print the collected operations.
        :return: List of dictionaries with 'type' and 'argument' keys.
        """
        operations = []

        # Simple greedy pattern - matches everything until the last ) on the line
        lattice_call_pattern = r"sim\.lattice\.(append|extend)\((.*)\)"
        
        # Store sim.lattice.append and sim.lattice.extend calls
        for match in re.finditer(lattice_call_pattern, content):
            operation, arg = match.groups()
            operations.append((match.start(), {"type": operation, "argument": arg.strip()}))

        # Store .reverse() calls
        reverse_pattern = r"(\w+)\.reverse\(\)"
        for match in re.finditer(reverse_pattern, content):
            operations.append((match.start(), {"type": "reverse", "argument": match.group(1)}))

        # important: sort operations by their position in the content
        # since the for loops can be executed in any order
        operations = [type for _, type in sorted(operations, key=lambda x: x[0])]

        if debug:
            self.print_lattice_operations(operations)

        return operations

    def _flatten(self, content: str, variable_name: str, debug: bool = True) -> List[str]:
        """
        Recursively expands a varible name to replace it with its set of elements.

        EX:
            content = '''
                cell = [drift1, quad1]
                line = [cell, cell]
                sim.lattice.extend(line)
            '''
            _flatten(content, "line")
                
            Results in:
                ["drift1", "quad1", "drift1", "quad1"]

        :param content: Full text content containing all variable definitions
        :param variable_name: Name of the specific variable to expand (e.g. "line")
        :param debug: Whether to print the expanded list.
        :return: List of individual element names with all nesting resolved.
        """

        # Check if the input is an inline list like "[monitor, elements.Drift(...)]"
        if variable_name.startswith("[") and variable_name.endswith("]"):
            list_contents = variable_name[1:-1]  # contents between brackets
            # split on commas that are NOT inside parentheses
            list_to_flatten = [element.strip() for element in re.split(r",\s*(?![^()]*\))", list_contents) if element.strip()]
        else:
            var_assignment_pattern = rf"{re.escape(variable_name)}\s*=\s*\[(.*?)\]"
            match = re.search(var_assignment_pattern, content, re.DOTALL)

            if not match:
                return [variable_name]  # It's not a list, it's a single element

            if debug:
                print(f"\nExpanding variable list definition for '{variable_name}':")
                print(f"  {match.group(0)}")

            list_content = match.group(1)
            list_to_flatten = [item.strip() for item in list_content.split(",") if item.strip()]

        expanded = []
        for item in list_to_flatten:
            # recursively expand each item
            sub_items = self._flatten(content, item, debug)
            expanded.extend(sub_items)

        return expanded


    def replace_variables(self, content: str, raw_lattice: List[str]) -> List[str]:
        """
        This function is called to simplify the lattice list by replacing variable names with their corresponding constructor calls.

        EX:
            (input)
                drift1 = elements.Drift(ds=1.0)
                quad1 = elements.Quad(k=0.5)
                raw_lattice = ["drift1", "quad1"]
            (output)    
                raw_lattice = ["elements.Drift(ds=1.0)", "elements.Quad(k=0.5)"]

        :param content: Full text content of the ImpactX simulation file.
        :param raw_lattice: List of lattice element variable names or constructor calls, e.g. ["drift1", "quad1"].
        :return: List with variable names replaced by their corresponding constructor calls.
        """
        element_mapping = {}

        ellement_assignment_pattern = r"^\s*(\w+)\s*=\s*(elements\.\w+\(.*?\))"
        all_element_assignments = re.findall(
            ellement_assignment_pattern, content, re.MULTILINE | re.DOTALL
        )

        for var_name, element in all_element_assignments:
            element_mapping[var_name] = element

        if not element_mapping:
            return raw_lattice

        # Replace each item in the list
        # later can be optimized by not iterating over the whole raw_lattice list
        return [element_mapping.get(item, item) for item in raw_lattice]
        
    def extract_used_variables(self, parsed_lattice: Dict[str, Any]) -> Set[str]:
        """
        Extracts used lattice variables from parsed lattice data.
        """
        used_variables = set()

        for element in parsed_lattice["lattice_elements"]:
            for value in element["parameters"].values():
                used_variables.add(value)

        return used_variables

    #-----------------------------------------------------------------------------
    # Debug methods
    # -----------------------------------------------------------------------------

    def print_lattice_operations(self, operations: List[Dict[str, str]]) -> None:
        """
        Prints all lattice operations in the order they appear.
        """

        print("Full lattice operation sequence (in order):")
        for operation in operations:
            print(f"  {operation['type']}({operation['argument']})")
