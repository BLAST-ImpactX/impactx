"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import plotly.graph_objects as go
from . import LatticeVisualizerElements as DrawElements
from . import LatticeVisualizerUtils as Utils
from ..... import setup_server

server, state, ctrl = setup_server()
draw = DrawElements()
utils = Utils()

def _error_plot(fig: go.Figure) -> go.Figure:
    """
    The error plot to display when there is an issue with the lattice data.

    :return: A plotly figure with an error message.
    """
    fig.add_annotation(
        text="Error: No lattice elements found or invalid data.",
        showarrow=False,
        font=dict(size=16, color="red"),
        align="center",
        xref="paper", yref="paper",
        x=0.5, y=0.5
    )
    fig.update_layout(
        title="Lattice Visualization - Error",
        xaxis=dict(title="X (m)", scaleanchor="y", scaleratio=1),
        yaxis=dict(title="Y (m)"),
        plot_bgcolor="white",
        margin=dict(l=30, r=30, t=40, b=40),
    )
    return fig


def lattice_visualizer():
    """
    Displays the lattice visualization using a plotly figure.
    Called every time the lattice list is modified.
    
    The current parameters which affect the visualization are:
    - `ds`: Length of the element (default 1.0 m)
    - `dx`: X offset of the element (default 0.0 m)
    - `dy`: Y offset of the element (default 0.0 m)
    - `rotation`: Rotation of the element in degrees (default 0.0)
    - `k`: Quadrupole strength (default 0.0)
       - Defocusing quads are colored red, focusing quads are blue
    - `phi`: Bend angle in radians (default 0.0)

    Shapes:
    - `drift`: Straight line representing drift space
    - `quadrupole`: Rectangle representing quadrupole magnets
    - `bend`: Arc representing bending magnets
    - `monitor`: Small rectangle representing beam monitors

    TO-DO:
    - `rc`: Radius of curvature for bends (default 0.0)
    - `name`: Name of the element (used for labeling)
    - `color`: Color of the element (default based on type)
    - `psi`: Pole face angle in radians (default 0.0)
    - `g`: gap parameter in m
    - more parameters as needed
    """
    fig = go.Figure()
    
    if not state.selected_lattice_list:
        return _error_plot(fig)
    
    try:
        x, y, rotation = 0, 0, 0
        draw.reset_legend()
        
        # Check if we should show labels based on element count
        element_count = len(state.selected_lattice_list)
        draw.set_show_labels(element_count <= 20)

        for index, element in enumerate(state.selected_lattice_list, 1):
            element_name = element.get("name", "").lower()  # This gets the actual element name
            element_label = utils.get_element_name_param(element)  # This gets the name from parameters
            ds = utils.get_element_param(element, "ds", 1.0)
            dx = utils.get_element_param(element, "dx", 0.0)
            dy = utils.get_element_param(element, "dy", 0.0)
            element_rotation = utils.get_element_param(element, "rotation", 0.0)
            rotation_total = rotation + element_rotation

            # Classify and draw element based on name
            if "drift" in element_name:
                x, y, rotation = draw.drift(fig, x, y, ds, dx, dy, rotation_total, element_label, index, element_name)
            elif "quad" in element_name:
                k = utils.get_element_param(element, "k", 0.0)
                x, y, rotation = draw.quad(fig, x, y, k, ds, dx, dy, rotation_total, element_label, index, element_name)
            elif "bend" in element_name or "dipole" in element_name:
                if element_name.startswith("sbend"):
                    rc = utils.get_element_param(element, "rc", 0.0)
                    x, y, rotation = draw.sBend(fig, x, y, ds, dx, dy, rotation_total, rc, element_label, index, element_name)
                elif element_name.startswith("exactsbend"):
                    phi = utils.get_element_param(element, "phi", 0)
                    x, y, rotation = draw.exactSBend(fig, x, y, ds, dx, dy, rotation_total, phi, element_label, index, element_name)
                else:
                    # Default bend handling
                    x, y, rotation = draw.drift(fig, x, y, ds, dx, dy, rotation_total, element_label, index, element_name)
            elif "monitor" in element_name or "bpm" in element_name:
                x, y, rotation = draw.beam_monitor(fig, x, y, rotation_total, ds, element_label, index, element_name)
            else:
                # Default fallback to drift
                x, y, rotation = draw.drift(fig, x, y, ds, dx, dy, rotation_total, element_label, index, element_name)
                
    except ValueError:
        return _error_plot(fig)

    fig.update_layout(
        title="Lattice Visualization",
        xaxis=dict(title="X (m)", scaleanchor="y", scaleratio=1),
        yaxis=dict(title="Y (m)"),
        plot_bgcolor="white",
        margin=dict(l=30, r=30, t=40, b=40),
    )
    
    return fig