"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from contextlib import contextmanager

from ... import html, setup_server, vuetify
from ...Input.components.card import CardComponents
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

    @staticmethod
    def icon_button(**kwargs):
        """
        Creates an icon with default properties
        specifically for the simulation history panels.
        """
        return CardComponents.card_button(
            density="default",
            size="small",
            **kwargs,
        )

    @staticmethod
    @contextmanager
    def view_details_card(title: str = "", prepend_icon: str = None, **kwargs):
        """
        Creates a card component used in the 'View Details'
        dialog of the simulation history.

        """

        with vuetify.VCard(
            rounded="lg",
            elevation=2,
            classes="pa-4 flex-grow-1",
            style="min-width: 150px;",
            **kwargs,
        ):
            with html.Div(classes="d-flex align-center mb-2"):
                if prepend_icon:
                    vuetify.VIcon(
                        prepend_icon, size="small", color="primary", classes="mr-2"
                    )
                html.Div(title, classes="text-caption font-weight-medium")
            yield

    @staticmethod
    @contextmanager
    def dialog(
        title: str,
        prepend_icon: str = None,
        dialog_var: str = None,
        width: str = "500px",
        **kwargs,
    ):
        """
        A dialog layout with a preset toolbar, icon, and close button.

        :param title: Title of the dialog.
        :param prepend_icon: Optional MDI icon name.
        :param dialog_var: State variable to close the dialog (e.g., 'sim_rename_dialog').
        :param width: Width of the dialog.
        :param kwargs: Additional keyword args passed to VCard.
        """
        with vuetify.VDialog(
            v_model=(dialog_var, False),
            max_width=width,
        ):
            with vuetify.VCard(
                elevation=10, classes="rounded-lg", style="overflow: hidden", **kwargs
            ):
                with vuetify.VToolbar(color="primary", classes="px-4"):
                    if prepend_icon:
                        vuetify.VIcon(prepend_icon)
                    vuetify.VToolbarTitle(title)
                    vuetify.VSpacer()
                    if dialog_var:
                        vuetify.VBtn(icon="mdi-close", click=f"{dialog_var} = false")
                yield
