"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

from .... import ctrl, html, vuetify
from ....Input.components.card import CardComponents


def _selector_chip(label: str, **kwargs) -> vuetify.VChip:
    """
    Renders a VChip for the header selector with customizable props.
    """
    return vuetify.VChip(
        label,
        key=("item",),
        classes="align-self-start mb-1",
        size="default",
        color="primary",
        variant="flat",
        **kwargs,
    )


class HeaderSelectorComponents:
    def header_chip(self, section_name: str, chip_type: str = "selected"):
        """
        Renders the VChip's shown in the selector.
        The properties of the chip depend on the section name.
        """

        label = f"{{{{ {section_name} }}}}"
        if chip_type == "selected":
            return _selector_chip(
                label,
                closable=("( !unselectable_headers.includes(item) )",),
                click_close=(ctrl.remove_selected_header, f"[{section_name}]"),
            )
        else:  # available
            return _selector_chip(
                label,
                append_icon="mdi-plus",
                click=(ctrl.add_selected_header, f"[{section_name}]"),
                style="cursor: pointer;",
            )

    def search_field(self) -> vuetify.VTextField:
        """
        Creates a VTextField with default properties
        for the header selection search.
        """

        return vuetify.VTextField(
            v_model=("over_s_header_search",),
            placeholder="Search",
            prepend_inner_icon="mdi-magnify",
            variant="outlined",
            density="compact",
            hide_details=True,
            clearable=True,
            classes="text-body-2",
        )

    def reset(self) -> vuetify.VBtn:
        """
        Creates a reset button for the header selection.
        Resets the selected headers to the default state.
        """

        return CardComponents.card_button(
            icon_name="mdi-refresh",
            click=ctrl.reset_selected_headers,
            description="Reset",
        )

    def _status_icon(self, v_model_name: str) -> vuetify.VIcon:
        """Displays an icon indicating the status of the section."""

        return vuetify.VIcon(
            "mdi-checkbox-blank-circle",
            size="x-small",
            color=(f"{v_model_name}.length > 0 ? 'success' : 'error'",),
            classes="mr-2",
        )

    def _count_chip(self, v_model_name: str) -> vuetify.VChip:
        """Displays the count of items in the section."""

        return vuetify.VChip(
            f"{{{{ {v_model_name}.length }}}}",
            size="x-small",
            color="primary",
            variant="flat",
        )

    def empty_state(self, icon: str, message: str):
        """
        Creates an empty state display.
        """
        with html.Div(
            classes="d-flex flex-column align-center justify-center h-100 pa-4"
        ):
            vuetify.VIcon(icon, size="48", color="grey-lighten-2", classes="mb-2")
            html.Span(message, classes="text-caption text-grey-darken-1 text-center")

    def section_header(self, name: str, v_model_name: str) -> None:
        with html.Div(
            classes="d-flex align-center justify-space-between pa-3 bg-grey-lighten-4"
        ):
            with html.Div(classes="d-flex align-center"):
                self._status_icon(v_model_name)
                html.Span(name, classes="text-body-2 font-weight-bold")
            self._count_chip(v_model_name)

    def selected_headers_column(self):
        """
        Displays the section for selected headers.
        """

        with vuetify.VCol(cols=6, classes="d-flex flex-column h-100"):
            components.section_header(name="Selected", v_model_name="selected_headers")
            vuetify.VDivider()
            with vuetify.Template(v_if="selected_headers.length == 0"):
                with html.Div(classes="pa-2", style="visibility: hidden"):
                    components.search_field()
                vuetify.VDivider(style="visibility: hidden")
            with html.Div(classes="pa-2 flex-grow-1"):
                with vuetify.Template(
                    v_if="selected_headers.length > 0",
                    v_for="item in selected_headers",
                ):
                    components.header_chip("item", "selected")
                with vuetify.Template(v_else=True):
                    components.empty_state("mdi-inbox-outline", "No headers selected")

    def available_headers_column(self):
        """
        Displays the section for available headers.
        """

        with vuetify.VCol(cols=6, classes="d-flex flex-column h-100"):
            components.section_header(
                name="Available", v_model_name="available_headers_on_ui"
            )
            vuetify.VDivider()
            with html.Div(classes="pa-2"):
                components.search_field()
            vuetify.VDivider()
            with html.Div(
                classes="pa-2 flex-grow-1", style="min-height: 0; overflow-y: auto;"
            ):
                with vuetify.Template(
                    v_if="available_headers_on_ui?.length > 0",
                    v_for="item in available_headers_on_ui",
                ):
                    components.header_chip("item", "available")

                with vuetify.Template(v_if="no_results"):
                    components.empty_state("mdi-magnify", "No matches found")
                with vuetify.Template(v_if="no_available_items"):
                    components.empty_state("mdi-inbox-outline", "No headers available")


components = HeaderSelectorComponents()
