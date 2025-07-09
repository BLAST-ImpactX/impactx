"""
This file is part of ImpactX

Copyright 2024 ImpactX contributors
Authors: Parthib Roy, Axel Huebl
License: BSD-3-Clause-LBNL
"""

import plotly.graph_objects as go

from ... import state


def over_s_plot():
    """
    Generates a 1D line plot using Plotly based on selected headers and filtered data.
    """

    selected_y = state.selected_headers
    over_s_data = state.over_s_table_data

    x = [row["s"] for row in over_s_data] if over_s_data else []
    y_axis = selected_y

    figure_data = []
    if y_axis:
        for column in y_axis:
            y = [row[column] for row in over_s_data]
            trace = go.Scatter(
                x=x,
                y=y,
                mode="lines+markers",
                name=column,
                line=dict(width=2),
                marker=dict(size=8),
            )
            figure_data.append(trace)

    return go.Figure(
        data=figure_data,
        layout=go.Layout(
            title="Over-S Plot",
            xaxis=dict(title="s"),
            yaxis=dict(title="m"),
            margin=dict(l=20, r=20, t=25, b=30),
        ),
    )
