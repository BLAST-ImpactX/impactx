"""
In situ visualization: the "survey3d" preset, the lattice, the orbit, the envelope and the
beam in 3D, rendered by Ascent

This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""

import math
import numbers
import os

from .InSitu import _COMPONENT_MOMENTS, _bin_range, _my_rank, _Visualization


class _Survey3D(_Visualization):
    """Ascent renders the lattice, the reference orbit, the envelope and the beam in 3D

    The beam is placed at the reference particle, in the global frame of the survey,
    with its coordinates and the envelope magnified.
    """

    #: opening angle of the camera [deg]
    FOV = 30.0
    #: opening angle of the racecar camera [deg]
    RACECAR_FOV = 60.0
    #: racecar camera: behind the end of the bunch (3 rendered rms bunch lengths behind
    #: its center) by this many element sizes
    RACECAR_GAP = 1.0
    #: racecar camera: above and to the side of the orbit, in element sizes, inside the
    #: bore (at least 0.44 element sizes in radius)
    RACECAR_ABOVE = 0.3
    RACECAR_SIDE = 0.2
    #: racecar camera: looks down the orbit this many element sizes ahead of the bunch
    RACECAR_AHEAD = 4.0

    def __init__(
        self,
        lattice,
        ref=None,
        envelope=None,
        magnify="auto",
        element_size="auto",
        camera="overview",
        image_size=(1280, 720),
        palette="cern-lhc",
        opacity=0.3,
        color="momentum_t",
        radius=None,
        samples=None,
        bore="round",
        cutaway=False,
        lighting=None,
        supersample=2,
        camera_smoothing="auto",
    ):
        if lattice is None:
            raise ValueError("survey3d: the lattice is needed (lattice=sim.lattice).")
        if envelope is not None and envelope < 0.0:
            raise ValueError("survey3d: envelope must be >= 0 rms beam sizes.")
        if not (camera in ("overview", "beam", "racecar") or isinstance(camera, dict)):
            raise ValueError(
                "survey3d: camera must be 'overview', 'beam', 'racecar' or Ascent "
                "camera parameters."
            )
        if not (magnify == "auto" or _positive(magnify)):
            raise ValueError(
                "survey3d: magnify must be 'auto', a factor, or a pair of factors "
                "(transverse, longitudinal)."
            )
        if not (element_size == "auto" or _positive(element_size)):
            raise ValueError("survey3d: element_size must be 'auto' or > 0 m.")
        if lighting not in (None, "camera", "cinematic"):
            raise ValueError(
                "survey3d: lighting must be None, 'camera' or 'cinematic'."
            )
        if bore not in ("round", "square", None):
            raise ValueError("survey3d: bore must be 'round', 'square' or None.")
        if not (isinstance(supersample, int) and supersample >= 1):
            raise ValueError("survey3d: supersample must be an integer >= 1.")
        if not (
            isinstance(opacity, numbers.Real)
            and math.isfinite(opacity)
            and 0.0 <= opacity <= 1.0
        ):
            raise ValueError("survey3d: opacity must be in [0, 1].")
        if not (
            camera_smoothing == "auto"
            or (
                isinstance(camera_smoothing, numbers.Real)
                and math.isfinite(camera_smoothing)
                and camera_smoothing >= 0.0
            )
        ):
            raise ValueError("survey3d: camera_smoothing must be 'auto' or >= 0 m.")
        self.lattice = lattice
        self.ref = ref
        self.envelope = envelope
        self.magnify = magnify
        self.element_size = element_size
        self.camera = camera
        self.image_size = image_size
        self.palette = palette
        self.opacity = opacity
        self.color = color
        self.radius = radius
        self.samples = samples
        self.bore = bore
        self.cutaway = cutaway
        self.lighting = lighting
        self.supersample = supersample
        self.camera_smoothing = camera_smoothing
        self._geometry = None
        self._meta = None
        # meshes of the elements, reused from frame to frame, see Geometry
        self._mesh_cache = {}

    def actions(self, pc=None, image_dir=None, element=None):
        import numpy as np
        import pandas as pd

        from ..plot.Survey import _survey
        from ..plot.Survey3D import (
            SURVEY_FIELD,
            WIRES_FIELD,
            frames,
            look_from,
            orbit_position,
            overview_axis,
            refined_lattice,
            smoothed_axes,
        )
        from ..plot.Survey3DGeometry import Geometry

        lattice = self.lattice
        if pc is not None and self.ref is not None:
            raise ValueError(
                "survey3d: pass the beam or a reference particle, not both."
            )
        ref = pc.ref if pc is not None else self.ref
        if ref is None:
            raise ValueError("survey3d: the beam or a reference particle is needed.")
        n_sigma = self.envelope
        if n_sigma is None:
            n_sigma = 1.0 if pc is not None else 0.0
        if n_sigma > 0.0 and pc is None:
            raise ValueError("survey3d: the envelope needs the beam.")

        # moments of the beam: MPI-collective
        covariance = pc.covariance_matrix() if n_sigma > 0.0 else None
        rbc = pc.beam_moments() if pc is not None else None
        num_particles = pc.total_number_of_particles() if pc is not None else 0

        # the index of the element: an element that is in the lattice more than once
        # is the one at the position of the reference particle
        envelope, index = _survey(lattice, ref, covariance=covariance, element=element)
        # sizes: elements relative to the extent of the lattice
        coarse = envelope[["x", "y", "z"]].to_numpy(dtype=np.float64)
        diagonal = float(np.linalg.norm(coarse.max(axis=0) - coarse.min(axis=0)))
        if self.element_size == "auto":
            a = 0.02 * diagonal if diagonal > 0.0 else 0.1
        else:
            a = float(self.element_size)

        # opaque elements, closed or cut open, are surfaces, shaded here: with a bore,
        # drifts are beam pipes
        opaque = self.opacity == 1.0
        surfaces = self.cutaway or opaque
        # the geometry: elements in cells of at most a quarter of their size along s,
        # for smooth lighting
        refined = refined_lattice(
            lattice,
            envelope,
            max_length=0.25 * a,
            pipes=surfaces and self.bore is not None,
        )
        survey = (
            envelope if refined is lattice else _survey(refined, ref, element=index)[0]
        )

        # magnification of the beam and the envelope: at its largest, the envelope
        # fills most of the bore of the elements, without an envelope 3 rms beam sizes
        # do; the bunch is at most a few elements long
        bore = Geometry.BORE * a
        if self.magnify == "auto":
            if n_sigma > 0.0:
                size = n_sigma * max(
                    envelope["sigma_x"].max(), envelope["sigma_y"].max()
                )
                m_perp = 0.9 * bore / size if size > 0.0 else 1.0
            elif rbc is not None and max(rbc["sigma_x"], rbc["sigma_y"]) > 0.0:
                m_perp = bore / (3.0 * max(rbc["sigma_x"], rbc["sigma_y"]))
            else:
                m_perp = 1.0
            length = 3.0 * ref.beta * rbc["sigma_t"] if rbc is not None else 0.0
            m_long = min(m_perp, 2.0 * a / length) if length > 0.0 else m_perp
        elif np.ndim(self.magnify) == 0:
            m_perp = m_long = float(self.magnify)
        else:
            m_perp, m_long = (float(m) for m in self.magnify)

        envelope_rows = envelope if n_sigma > 0.0 else None
        geometry = Geometry(
            refined,
            survey,
            envelope_rows,
            n_sigma,
            m_perp,
            a,
            self.palette,
            self.bore,
            self.cutaway,
            opaque,
            cache=self._mesh_cache,
        )
        self._geometry = geometry

        # the frame of the reference particle
        p_ref = np.array([ref.x, ref.y, ref.z])
        _, t, n, v = frames(
            pd.DataFrame(
                {
                    "x": [ref.x],
                    "y": [ref.y],
                    "z": [ref.z],
                    "px": [ref.px],
                    "py": [ref.py],
                    "pz": [ref.pz],
                }
            )
        )
        t, n, v = t[0], n[0], v[0]

        # camera
        width, height = self.image_size
        up = np.array([0.0, 1.0, 0.0])
        sigma = (
            8.0 * a if self.camera_smoothing == "auto" else float(self.camera_smoothing)
        )
        inside = False  # is the camera inside the lattice?
        if self.camera == "overview":
            # from above the side: the beam moves from left to right
            u = overview_axis(survey)
            side = np.array([-u[1], 0.0, u[0]])
            elevation = np.radians(35.0)
            direction = np.cos(elevation) * side + np.sin(elevation) * up
            camera = look_from(
                geometry.points(), direction, up, width / height, self.FOV
            )
        elif self.camera == "beam":
            # close to the beam, from above its side and a little behind: shows
            # the horizontal and the longitudinal extent. The camera follows the
            # orbit averaged over a few elements, so that it turns steadily.
            tc, nc, vc = smoothed_axes(survey, ref.s, sigma)
            azimuth, elevation = np.radians(20.0), np.radians(55.0)
            horizontal = np.cos(azimuth) * nc - np.sin(azimuth) * tc
            direction = np.cos(elevation) * horizontal + np.sin(elevation) * vc
            corners = np.array(
                [
                    p_ref + cn * 2.0 * a * nc + cv * 2.0 * a * vc + ct * 4.0 * a * tc
                    for cn in (-1.0, 1.0)
                    for cv in (-1.0, 1.0)
                    for ct in (-1.0, 1.0)
                ]
            )
            camera = look_from(corners, direction, up, width / height, self.FOV)
            # corners of the lattice behind the camera: see add_anchor
            geometry.add_anchor(camera["position"])
        elif self.camera == "racecar":
            # on the orbit behind the bunch, a little above and to the side of it,
            # looking down the beam pipe ahead. The view follows the orbit between
            # the camera and the point it looks at, averaged over a few elements, so
            # that it turns steadily.
            # the tail of the bunch: 3 rms lengths behind its center, or its last
            # particle if that is farther, e.g., in the tails of a Gaussian bunch
            tail = 0.0
            if rbc is not None:
                tail = ref.beta * m_long * max(3.0 * rbc["sigma_t"], rbc["max_t"], 0.0)
            behind = tail + self.RACECAR_GAP * a
            ahead = self.RACECAR_AHEAD * a
            tc, nc, vc = smoothed_axes(survey, ref.s + 0.5 * (ahead - behind), sigma)
            position = (
                orbit_position(survey, ref.s - behind)
                + self.RACECAR_ABOVE * a * vc
                + self.RACECAR_SIDE * a * nc
            )
            camera = {
                "position": position.tolist(),
                "look_at": (position + (behind + ahead) * tc).tolist(),
                "up": vc.tolist(),
                "fov": self.RACECAR_FOV,
                # inside the lattice: the default clipping range starts at a tenth
                # of the extent of the scene
                "near_plane": 0.01 * a,
                "far_plane": 10.0 * max(diagonal, a),
            }
            # surfaces of opaque elements next to the camera are drawn close up
            geometry.add_anchor(position)
            inside = True
        else:
            camera = self.camera

        # size of a pixel at the point the camera looks at
        if "position" in camera and "look_at" in camera:
            distance = np.linalg.norm(
                np.asarray(camera["position"]) - np.asarray(camera["look_at"])
            )
            fov = np.radians(camera.get("fov", self.FOV))
            pixel = 2.0 * distance * np.tan(fov / 2.0) / height
        else:
            pixel = 1.0e-3 * max(diagonal, a)
        # orbit and envelope: wires of about a pixel. With the camera inside the
        # lattice, as thin close to the camera as far away: else they would cover the
        # view
        geometry.add_wires(
            0.6 * pixel,
            ring_spacing=a,
            camera=(camera["position"], distance) if inside else None,
        )
        # a light at the camera
        dark = self.lighting == "cinematic"
        if dark and "position" in camera and "look_at" in camera:
            # a soft sky light, a hard key light and a spot light on the beam
            geometry.light_cinematic(camera["position"], camera["look_at"], focus=p_ref)
        elif (self.lighting == "camera" or surfaces) and "position" in camera:
            # opaque elements are always shaded: the renderer does not, see below
            geometry.light(camera["position"])

        plots = {}
        if geometry.num_colors > 0 and surfaces:
            # opaque elements: shaded by the ray tracer
            plots["lattice"] = {
                "type": "pseudocolor",
                "field": SURVEY_FIELD,
                "min_value": geometry.surface_value_range()[0],
                "max_value": geometry.surface_value_range()[1],
                "color_table": geometry.surface_color_table(),
            }
        elif geometry.num_colors > 0:
            plots["lattice"] = {
                "type": "volume",
                "field": SURVEY_FIELD,
                "min_value": 0.0,
                "max_value": geometry.max_value,
                "samples": self.samples
                if self.samples is not None
                else geometry.samples(),
                "color_table": geometry.color_table(self.opacity),
            }
        plots["wires"] = {
            "type": "pseudocolor",
            "field": WIRES_FIELD,
            "min_value": 0.0,
            "max_value": 1.0,
            "color_table": geometry.wires_color_table(
                orbit=(0.85, 0.85, 0.85) if dark else (0.0, 0.0, 0.0)
            ),
        }

        # the beam at the reference particle, magnified:
        # lab = p_ref + (n, v, t) diag(m_perp, m_perp, m_long) (x, y, -beta t)
        actions = []
        colorbar = None
        if num_particles > 0:
            col_long = -ref.beta * m_long * t
            matrix = [
                n[0] * m_perp, v[0] * m_perp, col_long[0], p_ref[0],
                n[1] * m_perp, v[1] * m_perp, col_long[1], p_ref[1],
                n[2] * m_perp, v[2] * m_perp, col_long[2], p_ref[2],
                0.0, 0.0, 0.0, 1.0,
            ]  # fmt: skip
            actions.append(
                {
                    "action": "add_pipelines",
                    "pipelines": {
                        "survey3d_placement": {
                            "f1": {
                                "type": "transform",
                                "params": {"topology": "beam", "matrix": matrix},
                            }
                        }
                    },
                }
            )
            radius = self.radius
            if radius is None:
                radius = max(0.01 * a, 0.6 * pixel)
            beam = {
                "type": "pseudocolor",
                "pipeline": "survey3d_placement",
                "field": self.color,
                "color_table": {"name": "Cool to Warm"},
                "points": {"radius": radius},
            }
            if self.color in _COMPONENT_MOMENTS:
                key, label = _COMPONENT_MOMENTS[self.color]
                vmin, vmax = _bin_range(
                    rbc[f"mean_{key}"] - 3.0 * rbc[f"sigma_{key}"],
                    rbc[f"mean_{key}"] + 3.0 * rbc[f"sigma_{key}"],
                )
                beam["min_value"] = vmin
                beam["max_value"] = vmax
                colorbar = (label, vmin, vmax)
            plots["beam"] = beam

        name = "survey3d" if image_dir is None else os.path.join(image_dir, "survey3d")
        render = {
            "image_name": name,
            # anti-aliasing: rendered larger, averaged down in frame()
            "image_width": width * self.supersample,
            "image_height": height * self.supersample,
            "annotations": "false",
            "bg_color": [0.06, 0.07, 0.09] if dark else [1.0, 1.0, 1.0],
            "fg_color": [0.9, 0.9, 0.9] if dark else [0.0, 0.0, 0.0],
            "camera": camera,
        }
        if surfaces:
            # the field values shade the elements smoothly: the renderer would shade
            # their flat faces
            render["shading"] = "disabled"
        actions.append(
            {
                "action": "add_scenes",
                "scenes": {"survey3d": {"plots": plots, "renders": {"r1": render}}},
            }
        )

        legend = list(geometry.legend.items())
        if geometry.has_envelope:
            from ..plot.Survey3D import ENVELOPE_COLOR

            legend.append((rf"envelope ({n_sigma:g}$\sigma$)", ENVELOPE_COLOR))
        note = None
        if num_particles > 0 or geometry.has_envelope:
            note = rf"beam magnified: transverse $\times${m_perp:.3g}" + (
                rf", longitudinal $\times${m_long:.3g}" if num_particles > 0 else ""
            )
        self._meta = {
            "legend": legend,
            "colorbar": colorbar,
            "note": note,
            "dark": dark,
        }
        return actions

    def geometry(self):
        if self._geometry is None or _my_rank() != 0:
            return None
        return self._geometry.blueprint()

    def frame(self, info, renderer, frame_name):
        frame = {"images": renderer.collect_images(info, frame_name)}
        if _my_rank() == 0:
            if self.supersample > 1:
                for path in frame["images"]:
                    _downsample(path, self.image_size)
            frame["survey3d"] = self._meta
        return frame


def _downsample(path, size):
    """Average an image file down to a size, in place: anti-aliasing by supersampling"""
    from PIL import Image

    try:
        box = Image.Resampling.BOX
    except AttributeError:  # Pillow < 9.1
        box = Image.BOX
    with Image.open(path) as im:
        small = im.convert("RGB").resize(tuple(size), box)
    small.save(path)


def _positive(value):
    """Is the value a positive number, or a sequence of positive numbers?"""
    import numpy as np

    try:
        values = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError):
        return False
    return (
        values.ndim <= 1
        and values.size in (1, 2)
        and bool(np.all(np.isfinite(values) & (values > 0.0)))
    )
