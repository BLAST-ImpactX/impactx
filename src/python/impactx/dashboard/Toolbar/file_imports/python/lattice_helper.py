import re

class DashboardLatticeConfigParser:

    @staticmethod
    def parse_lattice_elements(content: str) -> dict:
        """
        Parses lattice elements and also extracts used lattice variables.
        """
        lattice_order = DashboardLatticeConfigParser.collect_lattice_operations(content, debug=True)
        result = DashboardLatticeConfigParser.parse_cleaned_lattice(content)
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
        lattice_elements = re.findall(element_pattern, content)

        for element_name, parameter in lattice_elements:
            element = {
                "name": element_name,
                "parameters": {}
            }

            # Match parameters in the format key=value
            parameter_pattern = r"(\w+)=([^,\)]+)"  # EX: ds=1.0, k1=0.5
            all_parameters = re.findall(parameter_pattern, parameter)

            for parameter_name, value in all_parameters:
                element["parameters"][parameter_name] = value.strip()

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
