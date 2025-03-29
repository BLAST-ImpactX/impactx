"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""
from contextlib import contextmanager

from ...Input.components.input import InputComponents
from ... import html, setup_server, vuetify

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
    
    @contextmanager
    def sim_details_card(title: str = "", prepend_icon: str = None, **kwargs):
        """
        Creates a card component used in the 'View Details'
        dialog of the simulation history.

        """
        
        with vuetify.VCard(
            rounded="lg",
            elevation=2,
            classes="pa-4 flex-grow-1",
            style="min-width: 150px;",
            **kwargs
        ):
            with html.Div(classes="d-flex align-center mb-2"):
                if prepend_icon:
                    vuetify.VIcon(prepend_icon, size="small", color="primary", classes="mr-2")
                html.Div(title, classes="text-caption font-weight-medium")
            yield
            