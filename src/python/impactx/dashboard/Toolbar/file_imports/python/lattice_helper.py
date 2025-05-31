import re

class DashboardLatticeConfigParser:

    @staticmethod
    def parse_lattice_elements(content: str) -> dict:
        """
        Parses lattice elements and also extracts used lattice variables.
        """
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
    def extract_used_variables(parsed_lattice: dict) -> set:
        """
        Extracts used lattice variables from parsed lattice data.
        """
        used_variables = set()

        for element in parsed_lattice["lattice_elements"]:
            for value in element["parameters"].values():
                used_variables.add(value)

        return used_variables
