import re

class DashboardLatticeConfigParser:

    @staticmethod
    def parse_lattice_elements(content: str) -> dict:
        """
        Parses lattice elements and also extracts used lattice variables.
        """
        lattice_order = DashboardLatticeConfigParser.collect_lattice_operations(content, debug=True)

        expanded_elements = []
        for operation in lattice_order:
            operation_type = operation["type"]
            operation_arg = operation["argument"]

            if operation_type == "extend":
                expanded_elements.extend(DashboardLatticeConfigParser.flatten_variable_list_definition(content, operation_arg))

        clean_lattice_list = DashboardLatticeConfigParser.replace_variable_names_with_elements(content, expanded_elements)
        clean_lattice_list_str = '\n'.join(clean_lattice_list)
        result = DashboardLatticeConfigParser.parse_cleaned_lattice(clean_lattice_list_str)
        
        result["used_lattice_variables"] = DashboardLatticeConfigParser.extract_used_variables(result)

        return result

    @staticmethod
    def parse_cleaned_lattice(content: str) -> dict:
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


    @staticmethod
    def collect_lattice_operations(content: str, debug=False) -> list:
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

        lattice_call_pattern = r"sim\.lattice\.(append|extend)"
        pattern = re.compile(rf"{lattice_call_pattern}\((.*?)\)")


        # Store sim.lattice.append and sim.lattice.extend calls
        for match in pattern.finditer(content):
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
            DashboardLatticeConfigParser.print_lattice_operations(operations)

        return operations

    @staticmethod
    def flatten_variable_list_definition(content: str, variable_name: str, debug=True) -> list:
        """
        Recursively expands a varible name to replace it with its set of elements.

        EX:
            content = '''
                cell = [drift1, quad1]
                line = [cell, cell]
                sim.lattice.extend(line)
            '''
            flatten_variable_definition(content, "line")
                
            Results in:
                ["drift1", "quad1", "drift1", "quad1"]

        :param content: Full text content containing all variable definitions
        :param variable_name: Name of the specific variable to expand (e.g. "line")
        :param debug: Whether to print the expanded list.
        :return: List of individual element names with all nesting resolved.
        """


        var_assignment_pattern = rf"{re.escape(variable_name)}\s*=\s*\[(.*?)\]"
        match = re.search(var_assignment_pattern, content, re.DOTALL)

        if not match:
            return [variable_name]  # It's not a list, it's a single element

        list_content = match.group(1)
        items = [item.strip() for item in list_content.split(",") if item.strip()]

        expanded = []
        for item in items:
            # recursively expand each item
            sub_items = DashboardLatticeConfigParser.flatten_variable_list_definition(content, item, debug)
            expanded.extend(sub_items)

        if debug:
            print(f"\nExpanded list for {variable_name}:")
            print(f"  {expanded}")

        return expanded


    @staticmethod
    def replace_variable_names_with_elements(content: str, raw_lattice: list) -> list:
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
        
    @staticmethod
    def extract_used_variables(parsed_lattice: dict) -> set:
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

    @staticmethod
    def print_lattice_operations(operations: list) -> None:
        """
        Prints all lattice operations in the order they appear.
        """

        print("Full lattice operation sequence (in order):")
        for operation in operations:
            print(f"  {operation['type']}({operation['argument']})")
