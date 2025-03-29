"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from ... import setup_server, vuetify
from ...Input.components.input import InputComponents

server, state, ctrl = setup_server()

class SimulationHistoryComponents:

    @staticmethod
    def status_chip(obj_expr: str):
        """
        Renders a VChip for simulation status.
        """

        status_binding = f"{obj_expr}.status"
        
        return vuetify.VChip(
            f"{{{{ {status_binding} }}}}",
            color=(f"window.getSimStatusColor({status_binding})",),
            variant="elevated",
            size="small",
        )

    @staticmethod
    def text_field(**kwargs):
        """
        Creates a VTextField with default properties
        specifically for the simulation history panels.
        """
        
        return InputComponents.text_field(
            density="comfortable",
            hide_details=True,
            variant="outlined",
            input_type="text",
            **kwargs,
        )