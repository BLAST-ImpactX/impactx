"""
This file is part of ImpactX

Copyright 2025 ImpactX contributors
Authors: Parthib Roy
License: BSD-3-Clause-LBNL
"""

import plotly.graph_objects as go
import numpy as np

def transform(x, y, rotation_deg, dx):
    """
    Transform coordinates based on angle and displacement.
    """
    rotation_rad = np.radians(rotation_deg)
    x_new = x + dx * np.cos(rotation_rad)
    y_new = y + dx * np.sin(rotation_rad)
    return x_new, y_new

import numpy as np

def rotate_corners(x: float, y: float, rotation_deg: float, ds: float = 1.0, width: float = 0.1) -> np.ndarray:
    """
    Generates rectangle's corners after applying rotation matrix.
    This is utilized to properly visualize a rotated lattice element in Plotly.

    :param x: starting x-coordinate before the rotation
    :param y: starting y-coordinate before the rotation
    :param rotation_deg: Rotation angle in degrees, counterclockwise.
    :param ds: Length of the rectangle along the local X-axis (default is 1.0).
    :param width: Half of the rectangle's height (default is 0.1).
    :return: A NumPy array of shape (5, 2) with rotated (x, y) corner coordinates, closed for polygon plotting.

    """
    rotation_rad = np.radians(rotation_deg)

    corners = np.array([
        [0, -width],
        [ds, -width],
        [ds, width],
        [0, width],
        [0, -width]  # close polygon
    ])

    R = np.array([
        [np.cos(rotation_rad), -np.sin(rotation_rad)],
        [np.sin(rotation_rad),  np.cos(rotation_rad)],
    ])

    rotated = corners @ R.T + [x, y]
    return rotated
import numpy as np


class LatticeVisualizerElements:

    def __init__(self):
        self.seen_elements = set()
        self.show_labels = True  # Track whether to show labels

    def reset_legend(self):
        self.seen_elements.clear()

    def set_show_labels(self, show_labels):
        """Set whether to show element labels on the plot."""
        self.show_labels = show_labels

    def _add_to_legend(self, element_type):
        """
        Adds an element type to the legend if it hasn't been added already.
        """
        if element_type not in self.seen_elements:
            self.seen_elements.add(element_type)

    def _add_trace(self, fig, show_legend=False, legend_name=None, legend_group=None, **kwargs):
        """
        This is the function that actually draws on the plotly figure.
        """
        kwargs.setdefault("showlegend", show_legend)
        kwargs.setdefault("hoverinfo", "text")
        
        if show_legend and legend_name:
            kwargs.setdefault("name", legend_name)
            kwargs.setdefault("legendgroup", legend_group or legend_name)
        elif legend_group:
            kwargs.setdefault("legendgroup", legend_group)
            kwargs.setdefault("showlegend", False)
        
        # Properly handle hover information
        if 'text' in kwargs:
            kwargs['hovertemplate'] = kwargs['text'] + '<extra></extra>'
        else:
            kwargs.setdefault("hoverinfo", "skip")
            
        trace = go.Scatter(**kwargs)
        fig.add_trace(trace)

    def _generate_hover_text(self, label, ds, dx, dy, rotation, **special_params):
        """
        Generate standardized hover text for lattice elements.
        
        :param label: Element label/name
        :param ds: Element length
        :param dx: X displacement
        :param dy: Y displacement  
        :param rotation: Rotation angle
        :param special_params: Additional parameters specific to element type (k, rc, phi, etc.)
        :return: Formatted hover text string
        """
        text = (
            f"<b>{label}</b><br>"
            f"ds: {ds} m<br>"
        )
        
        # Add special parameters if provided
        for param_name, param_value in special_params.items():
            if param_name == 'k':
                text += f"k: {param_value} m⁻²<br>"
            elif param_name == 'rc':
                text += f"rc: {param_value} m<br>"
            elif param_name == 'phi':
                text += f"phi: {param_value}°<br>"
        
        text += (
            f"dx: {dx} m<br>"
            f"dy: {dy} m<br>"
            f"rotation: {rotation}°"
        )
        
        return text

    def drift(self, fig, x, y, ds, dx, dy, rotation, label):
        show_legend = self._add_to_legend("drift")
        rotation_rad = np.radians(rotation)
        thickness = 0.05  # line thickness (half-height for visual padding)
        x += dx
        y += dy

        rotated_corners = rotate_corners(x, y, rotation, ds, thickness)
        xs, ys = rotated_corners[:, 0], rotated_corners[:, 1]

        self._add_trace(
            fig,
            x=xs,
            y=ys,
            mode="lines",
            fill="toself",
            line=dict(color="gray", width=1),
            fillcolor="lightgray",
            text=self._generate_hover_text(label, ds, dx, dy, rotation),
            show_legend=show_legend,
            legend_name="Drift",
            legend_group="Drift"
        )

        if self.show_labels:
            self._add_annotation(
                fig,
                x=np.mean(xs),
                y=np.mean(ys) + 0.3,
                label=label,
                font=dict(size=10),
            )

        x1 = x + ds * np.cos(rotation_rad)
        y1 = y + ds * np.sin(rotation_rad)
        return x1, y1, rotation


    def quad(self, fig, x, y, k, ds, dx, dy, rotation, label):
        x1, y1 = transform(x, y, rotation, ds)
        x += dx
        y += dy
    
        match k:
            case _ if k > 0:
                quad_type = "focusing_quad"
                legend_name = "Focusing Quadrupole"
                line_color = "darkblue"
                fill_color = "lightblue"
            case _ if k < 0:
                quad_type = "defocusing_quad"
                legend_name = "Defocusing Quadrupole"
                line_color = "darkred"
                fill_color = "lightcoral"
            case _:
                quad_type = "quadrupole"
                legend_name = "Quadrupole"
                line_color = "darkgreen"
                fill_color = "lightgreen"

        show_legend = self._add_to_legend(quad_type)
        rotated_corners = rotate_corners(x, y, rotation, ds, 0.2)
        xs, ys = rotated_corners[:, 0], rotated_corners[:, 1]
    
        self._add_trace(
            fig,
            x=xs,
            y=ys,
            mode="lines",
            fill="toself",
            line=dict(color=line_color, width=2),
            fillcolor=fill_color,
            text=self._generate_hover_text(label, ds, dx, dy, rotation, k=k),
            show_legend=show_legend,
            legend_name=legend_name,
            legend_group=legend_name
        )

        if self.show_labels:
            self._add_annotation(fig, x=(x + x1)/2, y=y+0.4, label=label)
        return x1, y1, rotation

    def sBend(self, fig, x, y, ds, dx, dy, rotation, rc, label):
        show_legend = self._add_to_legend("sbend")
        """
        Draw a sector‐bend (SBEND) of length ds that has radius rc.
        - rc is the radius of curvature (in meters).
        - ds is the arc length (in meters).
        - rotation is the incoming reference angle (in degrees).
        - label is the magnet's name for hover/annotation.
        """

        # Apply any lateral offsets first:
        x += dx
        y += dy

        phi_rad = ds / rc           # Bend angle in radians
        rotation_rad = np.radians(rotation)

        # The circular‐arc center (in lab coords) is found by: 
        #   cx = x - ρ sin(incoming_angle)
        #   cy = y + ρ cos(incoming_angle)
        # because the bend is in the local x–z plane.
        r = rc
        cx = x - r * np.sin(rotation_rad)
        cy = y + r * np.cos(rotation_rad)

        n_points = 50
        arc_thetas = np.linspace(0, phi_rad, n_points)
        arc_x = cx + r * np.sin(rotation_rad + arc_thetas)
        arc_y = cy - r * np.cos(rotation_rad + arc_thetas)

        self._add_trace(
            fig,
            x=arc_x,
            y=arc_y,
            mode="lines",
            line=dict(color="blue", width=3),
            text=self._generate_hover_text(label, ds, dx, dy, rotation, rc=rc, phi=np.degrees(phi_rad)),
            show_legend=show_legend,
            legend_name="Sector Bend",
            legend_group="Sector Bend"
        )
        
        if self.show_labels:
            self._add_annotation(
                fig,
                x=np.mean(arc_x),
                y=np.mean(arc_y) + 0.3,
                label=label,
                font=dict(size=10),
            )

        # Compute exit point and exit angle (in degrees)
        x_end = arc_x[-1]
        y_end = arc_y[-1]
        final_angle = rotation + np.degrees(phi_rad)
        return x_end, y_end, final_angle

    def exactSBend(self, fig, x, y, ds: float, dx: float, dy: float, rotation_deg: float, phi_deg: float, label: str):
        show_legend = self._add_to_legend("exactsbend")
        """
        Draws an ExactSBend lattice element on the lattice visualization.
        """

        phi_rad = np.radians(phi_deg) # phi is given in degrees in the input
        rotation_rad = np.radians(rotation_deg) # may or may not be given

        #  Ensure the bend starts at the proper x and y coordinates
        x += dx
        y += dy

        # Compute curvature radius and arc path
        r = ds / phi_rad  # radius of curvature

        # find coords of center of the circle
        circle_center_x = x - r * np.sin(rotation_rad)
        circle_center_y = y + r * np.cos(rotation_rad)

        # Generate arc points
        n_points = 100 # determines smoothness of the arc
        arc_thetas = np.linspace(0, phi_rad, n_points)
        arc_x = circle_center_x + r * np.sin(rotation_rad + arc_thetas)
        arc_y = circle_center_y - r * np.cos(rotation_rad + arc_thetas)

        self._add_trace(
            fig,
            x=arc_x,
            y=arc_y,
            mode="lines",
            line=dict(color="blue", width=3),
            text=self._generate_hover_text(label, ds, dx, dy, rotation_deg, phi=phi_deg),
            show_legend=show_legend,
            legend_name="Exact Sector Bend",
            legend_group="Exact Sector Bend"
        )

        if self.show_labels:
            self._add_annotation(
                fig,
                x=np.mean(arc_x),
                y=np.mean(arc_y) + 0.3,
                label=label,
                font=dict(size=10),
            )

        # Compute new beamline exit point and angle
        final_angle = rotation_deg + phi_deg
        x_end = arc_x[-1]
        y_end = arc_y[-1]
        return x_end, y_end, final_angle

    def beam_monitor(self, fig, x, y, rotation, length, label):
        show_legend = self._add_to_legend("monitor")
        x1, y1 = transform(x, y, rotation, length)
        fig.add_shape(
            type="rect",
            x0=x, x1=x1,
            y0=y-0.15, y1=y+0.15,
            line=dict(color="darkgray"),
            fillcolor="lightgray",
        )
        self._add_trace(
            fig,
            x=[(x + x1)/2], y=[y],
            mode="markers",
            marker=dict(size=15, color='rgba(0,0,0,0)'),
            text=self._generate_hover_text(label, length, 0, 0, rotation),
            show_legend=show_legend,
            legend_name="Beam Monitor",
            legend_group="Beam Monitor"
        )
        if self.show_labels:
            self._add_annotation(fig, x=(x + x1)/2, y=y+0.3, label=label)
        return x1, y1, rotation

