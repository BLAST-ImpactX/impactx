"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

class LatticeVisualizerUtils:

    @staticmethod
    def get_element_param(element, name, default=0.0):
        for param in element.get("parameters", []):
            if param.get("parameter_name", "").lower() == name.lower():
                try:
                    return float(param.get("sim_input", default))
                except (ValueError, TypeError):
                    raise ValueError(f"Invalid value for {name}: {param.get('sim_input', default)}")
        return default

    @staticmethod
    def get_element_name_param(element):
        """Get the name parameter from the element's parameters list."""
        for param in element.get("parameters", []):
            if param.get("parameter_name", "").lower() == "name":
                name_value = param.get("sim_input", "")
                if name_value and name_value.strip():  # Check if name is not empty
                    return name_value
        
        # If no name parameter found or it's empty, create a meaningful fallback
        element_type = element.get("name", "")
        return f"{element_type}"  # Use the element type as fallback
