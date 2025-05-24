"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from .. import ctrl, html, state, vuetify
from ..Input.components.card import CardComponents
from ..Input.utils import GeneralFunctions
from ..Run.simulation import dashboard_sim_inputs

state.expand_all_sections = False


class InputToolbar:
    """
    Contains toolbar components for the 'Input' page.
    """

    @ctrl.trigger("export")
    def on_export_click():
        """
        Called when the export button is clicked.
        """
        return dashboard_sim_inputs(is_exporting=True)

    @ctrl.add("reset_all")
    def on_reset_all_click():
        """
        Called when the reset button is clicked.
        """
        GeneralFunctions.reset_inputs("all")

    @staticmethod
    def export_button() -> vuetify.VBtn:
        """
        Displays a button to export the current simulation inputs
        to a Python file.

        On click, it triggers the export of the current
        simulation inputs to a file named 'impactx_simulation.py'.
        """

        return vuetify.VBtn(
            "Export",
            click="utils.download('impactx_simulation.py', trigger('export'), 'text/plain')",
            variant="outlined",
            size="small",
            disabled=("disableRunSimulationButton", True),
            classes="mx-2",
            prepend_icon="mdi-download",
            color="#00313C",
        )

    @staticmethod
    def collapse_all_sections_button() -> vuetify.VBtn:
        """
        Displays a button to collapse or expand all input sections.

        On click, it toggles the visibility of all input sections
        in the dashboard, allowing users to quickly hide or show
        all input fields.
        """

        return CardComponents.card_button(
            ["mdi-collapse-all", "mdi-expand-all"],
            click=ctrl.collapse_all_sections,
            dynamic_condition="expand_all_sections",
            description=["Minimize All", "Show All"],
            classes="mr-2",
        )

    @staticmethod
    def import_button() -> vuetify.VBtn:
        """
        Displays a button to upload a file.

        On click, it opens a file input dialog to select
        an ImpactX simulation file (.py) for import.
        """

        vuetify.VFileInput(
            v_model=("import_file",),
            accept=".py",
            __properties=["accept"],
            style="display: none;",
            ref="fileInput",
        )
        with html.Div(style="position: relative;"):
            with vuetify.VBtn(
                "Import",
                click="trame.refs.fileInput.click()",
                size="small",
                variant="outlined",
                prepend_icon="mdi-upload",
                disabled=("(import_file_details)",),
                color=("import_file_error ? 'error' : '#00313C'",),
            ):
                pass
            with html.Div(
                style="position: absolute; font-size: 10px; width: 100%; padding-top: 2px; display: flex; justify-content: center; white-space: nowrap;"
            ):
                html.Span(
                    "{{ import_file_error ? import_file_error_message : import_file_details }}",
                    style="text-overflow: ellipsis; overflow: hidden;",
                    classes=(
                        "import_file_error ? 'error--text' : 'grey--text text--darken-1'",
                    ),
                )
                vuetify.VIcon(
                    "mdi-close",
                    x_small=True,
                    style="cursor: pointer;",
                    click=ctrl.reset_all,
                    v_if="import_file_details || import_file_error",
                    color=("import_file_error ? 'error' : 'grey darken-1'",),
                )

    @staticmethod
    def reset_inputs_button() -> vuetify.VBtn:
        """
        Displays a button to reset all input fields.

        On click, it resets all input fields to their
        default values.
        """

        return vuetify.VBtn(
            "Reset",
            id="reset_all_inputs_button",
            click=ctrl.reset_all,
            variant="outlined",
            size="small",
            prepend_icon="mdi-refresh",
            classes="mr-2",
            color="#00313C",
        )

    @staticmethod
    def error_notification():
        """
        Displays a notification when there are input validation errors.
        Shows error count and provides tooltip with details.
        """
        num_input_errors = "simulation_error_details.length"
        with html.Div(
            v_if="disableRunSimulationButton",
            style="display: flex; align-items: center; margin-right: 8px; padding: 4px 8px; background-color: #fff3e0; border-radius: 4px; border-left: 3px solid #ff9800;"
        ):
            html.Span(
                "{{ simulation_error_details.length }} input error{{ simulation_error_details.length > 1 ? 's' : '' }}",
                style="font-size: 12px; color: #e65100; margin-right: 4px;"
            )
            with vuetify.VTooltip(location="bottom"):
                with html.Template(v_slot_activator="{ props }"):
                    vuetify.VIcon(
                        "mdi-information",
                        size="small",
                        color="#ff9800",
                        v_bind="props",
                    )
                with html.Div():
                    with html.Div(
                        v_for="(category_group, category_index) in categorized_error_details",
                    ):
                        html.Div(
                            "{{ category_group.category }}",
                            style="color: #e65100;"
                        )
                        with html.Ul(style="padding-left: 16px;"):
                            with html.Li(
                                v_for="(error, error_index) in category_group.errors",
                                style="font-size: 12px;"
                            ):
                                html.Span("{{ error }}")