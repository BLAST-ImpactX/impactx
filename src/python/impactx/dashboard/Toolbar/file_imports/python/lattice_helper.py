import re

class DashboardLatticeConfigParser:

    @staticmethod
    def parse_lattice_elements(content: str) -> dict:
        """
        Parses lattice elements and also extracts used lattice variables.
        """
        result = DashboardLatticeConfigParser.parse_simple_lattice(content)
        result["used_lattice_variables"] = DashboardLatticeConfigParser.extract_used_variables(result)

        return result

    @staticmethod
    def parse_simple_lattice(content: str) -> dict:
        """
        Parses simple lattice elements from the simulation file content.
        """

        dictionary = {"lattice_elements": []}
        lattice_elements = re.findall(r"elements\.(\w+)\((.*?)\)", content)

        for element_name, element_parameter in lattice_elements:
            element = {"element": element_name, "parameters": {}}

            parameter_pairs = re.findall(r"(\w+)=([^,\)]+)", element_parameter)
            for parameter_name, parameter_value in parameter_pairs:
                parameter_value_cleaned = parameter_value.strip("'\"")
                element["parameters"][parameter_name] = parameter_value_cleaned

            dictionary["lattice_elements"].append(element)

        return dictionary

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
