
import re

class DashboardLatticeConfigParser:

    @staticmethod
    def parse_lattice_elements(content: str) -> dict:
        """
        Parses lattice elements from the simulation file content.

        :param content: The content of the ImpactX simulation file.
        """

        dictionary = {"lattice_elements": []}
        used_variables = set()

        lattice_elements = re.findall(r"elements\.(\w+)\((.*?)\)", content)

        for element_name, element_parameter in lattice_elements:
            element = {"element": element_name, "parameters": {}}

            parameter_pairs = re.findall(r"(\w+)=([^,\)]+)", element_parameter)
            for parameter_name, parameter_value in parameter_pairs:
                parameter_value_cleaned = parameter_value.strip("'\"")
                element["parameters"][parameter_name] = parameter_value_cleaned
                used_variables.add(parameter_value_cleaned)

            dictionary["lattice_elements"].append(element)

        dictionary["used_lattice_variables"] = used_variables
        return dictionary