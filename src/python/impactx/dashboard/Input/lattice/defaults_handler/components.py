"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

from .... import ctrl, html, vuetify


def header_summary():
    """
    Renders a small legend for match count next to search.
    """
    vuetify.VIcon(
        "mdi-checkbox-blank-circle", size="x-small", color="primary", classes="mr-2"
    )
    html.Span(
        "Matches: {{ lattice_defaults_filtered.length }} / {{ lattice_defaults.length }}",
        classes="text-caption text-grey-darken-1",
        aria_live="polite",
        __properties=["aria-live"],
    )


def text_field(**kwargs):
    """
    Shared text field with common defaults for this view.
    Distinct props can be passed via kwargs and will override defaults.
    """
    defaults = dict(
        variant="outlined",
        density="compact",
        hide_details=True,
        style="min-width: 0;",
    )
    defaults.update(kwargs)
    return vuetify.VTextField(**defaults)


def search_bar():
    """
    Search row with text field and match summary.
    """
    with vuetify.VRow(classes="align-center"):
        with vuetify.VCol(cols=True):
            text_field(
                v_model=("lattice_defaults_filter", ""),
                label="Search parameters",
                placeholder="e.g., nslice",
                prepend_inner_icon="mdi-magnify",
                clearable=True,
                classes="text-body-2",
                id="lattice_defaults_search",
                aria_label="Search parameters",
                __properties=["aria-label"],
            )
        with vuetify.VCol(cols="auto", classes="d-flex align-center ml-3"):
            header_summary()


def pagination():
    """
    Pagination control bound to lattice defaults filter list.
    """
    return vuetify.VPagination(
        v_model=("lattice_defaults_page", 1),
        length=("Math.max(1, Math.ceil(lattice_defaults_filtered.length / 5))",),
        total_visible=7,
        __properties=["length", "total_visible"],
        density="comfortable",
    )


def reset_button():
    """
    Reset defaults button.
    """
    return vuetify.VBtn(
        "Reset Defaults",
        id="reset_lattice_defaults",
        color="primary",
        click=ctrl.reset_lattice_defaults,
        block=True,
    )
