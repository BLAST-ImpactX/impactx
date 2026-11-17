"""
The lattice, the envelope and the orbit of the 3D survey as a Conduit Mesh Blueprint

This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""

import math

from .Survey3D import (
    _CAGE_WIRES,
    _RING_VERTICES,
    ENVELOPE_COLOR,
    PIPE_KINDS,
    SKIP_KINDS,
    SURVEY_FIELD,
    SURVEY_TOPOLOGY,
    WIRE_ENVELOPE,
    WIRE_ORBIT,
    WIRES_FIELD,
    WIRES_TOPOLOGY,
    frames,
)
from .Survey3DMeshes import (
    _add_hexahedra,
    _box_sweep_varying,
    _normalized,
    _row_dots,
    _row_norms,
    _shade,
    _shell_sweep,
    _smoothstep,
    _surface_patches,
    _unchanged,
    _unique_positions,
)


def _cross_section(kind):
    """Half width and half height of an element, relative to the element size"""
    if "bend" in kind:
        return 1.5, 0.7
    if "Quad" in kind:
        return 1.0, 1.0
    if "RF" in kind or "Buncher" in kind or "Cavity" in kind:
        return 1.2, 1.2
    return 0.9, 0.9


class Geometry:
    """The lattice, the envelope and the orbit as one Conduit Mesh Blueprint domain

    The elements are boxes, by default with a round bore that holds the envelope, as a magnet
    holds the beam pipe. They are one hexahedral topology for one semi-transparent
    volume plot. The field ``survey_kind`` encodes the color at each vertex: for the
    color index k = 1 ... number of colors, the walls have values in [k, k + SPAN]
    and the lining of the bore in [k + 0.5, k + 0.5 + SPAN], from unlit to lit, see
    light(). Opaque elements, closed or cut open, are their surfaces instead, one
    quadrilateral topology, with values in [k, k + SPAN]; with a bore, drifts are beam
    pipes around the same bore.

    The orbit and the envelope are thin opaque wires in a second hexahedral
    topology, visible through the elements: the envelope is a cage of longitudinal
    wires and rings. The field ``survey_wire_kind`` is WIRE_ORBIT or WIRE_ENVELOPE.

    Parameters
    ----------
    lattice : KnownElementsList
    geometry : pandas.DataFrame
        Survey of the lattice for the geometry, one row per (refined) slice.
    envelope : pandas.DataFrame, optional
        Survey with the rms beam sizes, one row per slice of the lattice.
    n_sigma : float
        Envelope half size in rms beam sizes.
    magnify : float
        Transverse magnification of the envelope.
    element_size : float
        Half size of the elements [m].
    palette : str
        Color palette of the elements, see ElementColors.
    bore : str, optional
        ``"round"`` or ``"square"`` bore through the elements, or None for solid ones.
    cutaway : bool
        Opaque elements, only their lower half, cut open through the beam axis.
    opaque : bool
        Opaque elements, closed around the bore.
    cache : dict, optional
        Meshes of the elements of an earlier geometry, e.g., of the previous frame of an
        animation: an element whose position, axes, size and bore did not change reuses
        its mesh. Updated in place; empty at first.
    """

    #: half size of the bore of the elements, relative to the element size
    BORE = 0.5
    #: least wall thickness of the elements, relative to the element size
    WALL = 0.4
    #: thickness of zero-length elements, relative to the element size
    THIN = 0.3
    #: opaque elements: wall thickness of the beam pipes in drifts, relative to the
    #: element size, their color and their legend label
    PIPE_WALL = 0.08
    PIPE_COLOR = (0.6, 0.61, 0.64)
    PIPE_LABEL = "beam pipe"
    #: opaque elements: zero-length elements reach this far into the bore, relative to
    #: the element size: a collar, that does not share surfaces with its neighbors
    COLLAR = 0.04
    #: an element reuses its mesh from the cache if its position and size changed by at
    #: most this much, relative to the element size, and its axes by at most this much
    CACHE_TOLERANCE = 1.0e-9
    #: thickness of the darker lining of the bore, relative to the element size
    LINER = 0.06
    #: brightness of the lining of the bore, relative to the element color
    LINER_SHADE = 0.4
    #: opacity of the lining of the bore, per crossing: rays cross it twice
    LINER_OPACITY = 0.2
    #: thickness of the layer between the wall and the lining of the bore, relative to the
    #: element size: wall and lining have their own vertices, for their own shading
    SEAM = 0.01
    #: radius of the bore through solid translucent elements, relative to the element
    #: size: too small to see, it gives the cells radial normals for the lighting
    PINHOLE = 0.02
    #: brightness of the lining of the bore with lighting, relative to the element color
    LINER_SHADE_LIT = 0.8
    #: brightness of surfaces that face away from the light
    AMBIENT = 0.5
    #: brightness of the surfaces of opaque elements that face away from the light
    OPAQUE_AMBIENT = 0.35
    #: light both sides of a surface: the back walls of translucent elements, too
    TWO_SIDED = True
    #: range of field values per color, from unlit to lit
    SPAN = 0.45
    #: cinematic lighting: weights of the soft sky light, the hard key light, the fill
    #: light from below and the spot light on the beam
    SKY, KEY, FILL, SPOT = 0.22, 0.28, 0.12, 0.8
    #: cinematic lighting: color of the shadows and of the highlights, added to the
    #: element color at the dark and at the bright end
    SHADOW_TINT = (0.02, 0.03, 0.07)
    HIGHLIGHT_TINT = (0.15, 0.13, 0.10)

    def __init__(
        self,
        lattice,
        geometry,
        envelope,
        n_sigma,
        magnify,
        element_size,
        palette,
        bore="round",
        cutaway=False,
        opaque=False,
        cache=None,
    ):
        import numpy as np
        from matplotlib.colors import to_rgb

        from ..extensions.KnownElementsList import _element_kind
        from .ElementColors import get_element_color

        if bore not in ("round", "square", None):
            raise ValueError("survey 3D: bore must be 'round', 'square' or None.")
        self.element_size = element_size
        p, t, n, v = frames(geometry)
        element = geometry["element"].to_numpy()
        s = geometry["s"].to_numpy(dtype=np.float64)

        # one frame per position: the envelope and the orbit
        rings = _unique_positions(s, np.arange(len(s)))
        self._p, self._t, self._n, self._v = p[rings], t[rings], n[rings], v[rings]
        self._s = s[rings]

        # the envelope: the beam sizes between the slice boundaries are interpolated
        a = element_size
        self._rx = self._ry = None
        if envelope is not None:
            s_env, keep = np.unique(envelope["s"].to_numpy()[::-1], return_index=True)
            keep = len(envelope) - 1 - keep  # the last row at each s
            sx = envelope["sigma_x"].to_numpy()[keep]
            sy = envelope["sigma_y"].to_numpy()[keep]
            self._rx = n_sigma * magnify * np.interp(self._s, s_env, sx)
            self._ry = n_sigma * magnify * np.interp(self._s, s_env, sy)
            r_tube = np.maximum(self._rx, self._ry)
        else:
            r_tube = np.zeros(len(self._s))

        def largest_tube(s0, s1):
            """Largest half size of the envelope in [s0, s1]"""
            # the positions are in increasing order
            first = np.searchsorted(self._s, s0, side="left")
            last = np.searchsorted(self._s, s1, side="right")
            ends = np.interp([s0, s1], self._s, r_tube)
            return float(max(np.max(r_tube[first:last], initial=0.0), *ends))

        coords = []  # arrays of (vertices, 3)
        conn = []  # arrays of (cells, 8)
        kinds = []  # color index per vertex
        linings = []  # per vertex: does it line the bore?
        normals = []  # per vertex: outward normal of its surface, zero if none
        num_vertices = 0

        #: legend entries: label -> RGB color
        self.legend = {}
        self.colors = []  # unique RGB colors, index + 1 is the field value
        self._walls = {}  # color index -> thinnest wall [m]
        self._liner = None  # thickness of the lining of the bores [m]

        indices = {}  # label -> color index

        def color_index(kind, rgb=None):
            if kind in indices:
                return indices[kind]
            if rgb is None:
                rgb = tuple(to_rgb(get_element_color(kind, palette=palette)))
            if rgb not in self.colors:
                self.colors.append(rgb)
            self.legend.setdefault(kind, rgb)
            indices[kind] = self.colors.index(rgb) + 1
            return indices[kind]

        # opaque elements, closed or cut open: their surfaces
        surfaces = cutaway or opaque
        if cache is not None:
            # meshes of other settings do not apply; the element size within roundoff,
            # e.g., a size relative to the extent of the lattice
            settings = (
                bore,
                cutaway,
                opaque,
                self.BORE,
                self.WALL,
                self.THIN,
                self.PIPE_WALL,
                self.COLLAR,
                self.LINER,
                self.SEAM,
                self.PINHOLE,
            )
            size = cache.get("element_size")
            if (
                size is None
                or abs(size - a) > self.CACHE_TOLERANCE * a
                or cache.get("settings") != settings
            ):
                cache.clear()
                cache["element_size"] = a
                cache["settings"] = settings
        # the rows of each element: a survey lists them in the order of the lattice
        if np.any(np.diff(element) < 0):
            raise ValueError("survey 3D: the survey rows are not in lattice order.")
        first_rows = np.searchsorted(element, np.arange(len(lattice)), side="left")
        last_rows = np.searchsorted(element, np.arange(len(lattice)), side="right")
        for index, el in enumerate(lattice):
            kind = _element_kind(el)
            # opaque elements with a bore: drifts are beam pipes
            pipe = surfaces and bore is not None and kind in PIPE_KINDS
            if kind in SKIP_KINDS and not pipe:
                continue
            # each read of a property is a call into C++
            zero_length = el.ds == 0.0
            pipe = pipe and not zero_length
            if kind in SKIP_KINDS and not pipe:
                continue  # a zero-length drift
            rows = np.arange(first_rows[index], last_rows[index])
            if not zero_length:
                # the row before the first slice is the entrance of the element
                rows = _unique_positions(s, np.concatenate(([rows[0] - 1], rows)))
                ends, tt, nn, vv = p[rows], t[rows], n[rows], v[rows]
                s0, s1 = s[rows[0]], s[rows[-1]]
            else:
                # a zero-length element: a thin slab
                r = rows[-1]
                half = 0.5 * self.THIN * a
                ends = np.array([p[r] - half * t[r], p[r] + half * t[r]])
                tt = np.array([t[r], t[r]])
                nn = np.array([n[r], n[r]])
                vv = np.array([v[r], v[r]])
                s0 = s1 = s[r]
            # the envelope fits into the bore, with walls around it
            radius = max(self.BORE * a, 1.1 * largest_tube(s0, s1))
            hw, hh = (f * a for f in _cross_section(kind))
            hw = max(hw, radius + self.WALL * a)
            hh = max(hh, radius + self.WALL * a)
            # the mesh: reused if position, axes, size and bore did not change
            lengths = np.concatenate((ends.ravel(), [radius, hw, hh]))
            axes = np.concatenate((tt.ravel(), nn.ravel(), vv.ravel()))
            cached = cache.get(index) if cache is not None else None
            if (
                cached is not None
                and cached[0] == (kind, pipe, zero_length)
                and _unchanged(cached[1], lengths, self.CACHE_TOLERANCE * a)
                and _unchanged(cached[2], axes, self.CACHE_TOLERANCE)
            ):
                mesh = cached[3]
            else:
                mesh = self._element_mesh(
                    ends,
                    tt,
                    nn,
                    vv,
                    radius,
                    hw,
                    hh,
                    bore,
                    pipe=pipe,
                    surfaces=surfaces,
                    cutaway=cutaway,
                    zero_length=zero_length,
                )
                if cache is not None:
                    cache[index] = ((kind, pipe, zero_length), lengths, axes, mesh)
            points, cells, lined, vertex_normals, wall, liner = mesh
            if liner > 0.0:
                self._liner = liner
            value = (
                color_index(self.PIPE_LABEL, self.PIPE_COLOR)
                if pipe
                else color_index(kind)
            )
            self._walls[value] = min(self._walls.get(value, wall), wall)
            coords.append(points)
            conn.append(cells + num_vertices)
            # per vertex: the shading is interpolated over the vertices of each cell
            kinds.append(np.full(len(points), value))
            linings.append(lined)
            normals.append(vertex_normals)
            num_vertices += len(points)

        self.num_colors = len(self.colors)
        self.max_value = float(max(self.num_colors, 1)) + 0.5 + self.SPAN
        self._surfaces = surfaces
        self._coords = np.concatenate(coords) if coords else np.zeros((0, 3))
        self._conn = np.concatenate(conn) if conn else np.zeros((0, 8), dtype=int)
        self._kinds = np.concatenate(kinds) if kinds else np.zeros(0, dtype=int)
        self._lined = np.concatenate(linings) if linings else np.zeros(0, dtype=bool)
        self._normals = np.concatenate(normals) if normals else np.zeros((0, 3))
        self._lit = False
        self._cinematic = False
        self.light(None)
        self._wires = None
        self._anchor = None

    def _element_mesh(
        self,
        ends,
        tt,
        nn,
        vv,
        radius,
        hw,
        hh,
        bore,
        pipe,
        surfaces,
        cutaway,
        zero_length,
    ):
        """The mesh of an element: points, cells, lining and normal per point, the
        thinnest wall and the thickness of the lining (0 for none)"""
        import numpy as np

        a = self.element_size
        liner = 0.0
        if pipe:
            # a round pipe, around the same bore as the elements
            outside = radius + self.PIPE_WALL * a
            points, cells, vertex_normals = _surface_patches(
                ends,
                tt,
                nn,
                vv,
                outside,
                outside,
                radius,
                "round",
                outside="round",
                lower_half=cutaway,
            )
            wall = self.PIPE_WALL * a
            lined = np.zeros(len(points), dtype=bool)
        elif surfaces:
            # opaque: the surfaces, as quadrilaterals
            if bore is None:
                radius = 0.0
            elif zero_length:
                radius -= self.COLLAR * a
            points, cells, vertex_normals = _surface_patches(
                ends,
                tt,
                nn,
                vv,
                hw,
                hh,
                radius,
                bore,
                lower_half=cutaway,
            )
            wall = min(hw, hh)
            lined = np.zeros(len(points), dtype=bool)
        else:
            # translucent elements are lined. Solid ones have a pinhole along the axis,
            # for radial normals.
            solid = bore is None
            liner = 0.0 if solid else self.LINER * a
            points, cells, lined, vertex_normals = _shell_sweep(
                ends,
                nn,
                vv,
                hw,
                hh,
                self.PINHOLE * a if solid else radius,
                "round" if solid else bore,
                liner,
                seam=self.SEAM * a,
            )
            # solid: as opaque through its center as the two walls around a bore
            wall = min(hw, hh) if solid else min(hw, hh) - radius - liner
        return points, cells, lined, vertex_normals, wall, liner

    def light(self, position):
        """Shade the elements for a light at a position, e.g., the camera

        Lambert shading: a surface has the brightness AMBIENT + (1 - AMBIENT) times the
        cosine of the angle between its normal and the direction to the light. With
        TWO_SIDED, the absolute cosine: the back walls of translucent elements are lit
        as well, else they have the brightness AMBIENT. The volume renderer has no lighting: the brightness
        is encoded in the field values, see color_table.

        Parameters
        ----------
        position : array_like or None
            Position of the light [m]. None: no lighting, all surfaces are lit.
        """
        import numpy as np

        if self._surfaces:
            # opaque elements, one-sided: without the flat shading of the renderer
            self._cinematic = False
            self._lit = position is not None
            lit = np.ones(len(self._kinds))
            if self._lit and len(self._kinds):
                to_light = _normalized(
                    np.asarray(position, dtype=float)[None, :] - self._coords
                )
                lit = np.clip(_row_dots(self._normals, to_light), 0.0, 1.0)
            self._values = self._kinds + self.SPAN * lit
            return
        self._lit = position is not None
        self._cinematic = False
        lit = np.ones(len(self._kinds))
        if self._lit and len(self._kinds):
            centers = self._coords  # lit at the vertices
            to_light = _normalized(np.asarray(position, dtype=float)[None, :] - centers)
            cosine = _row_dots(self._normals, to_light)
            if self.TWO_SIDED:
                cosine = np.abs(cosine)
            has_normal = _row_dots(self._normals, self._normals) > 0.0
            lit = np.where(has_normal, np.clip(cosine, 0.0, 1.0), 1.0)
        self._values = self._kinds + np.where(self._lined, 0.5, 0.0) + self.SPAN * lit

    def light_cinematic(self, camera_position, look_at, focus):
        """Shade the elements as in a film set, for a dark background

        A soft sky light (brighter from above), a hard key light from the upper left in
        front, and a spot light on the focus, e.g., the beam; distant elements fade
        into the background. The brightness is encoded in the field values, see
        color_table.

        Parameters
        ----------
        camera_position, look_at : array_like
            The camera [m].
        focus : array_like
            Where the spot light points to [m].
        """
        import numpy as np

        self._lit = True
        self._cinematic = True
        if not len(self._kinds):
            self._values = np.zeros(0)
            return
        a = self.element_size
        up = np.array([0.0, 1.0, 0.0])
        camera = np.asarray(camera_position, dtype=float)
        look_at = np.asarray(look_at, dtype=float)
        focus = np.asarray(focus, dtype=float)
        back = _normalized((camera - look_at)[None, :])[0]  # towards the camera
        right = _normalized(np.cross(up, back)[None, :])[0]

        normals = self._normals
        has_normal = _row_dots(normals, normals) > 0.0
        centers = self._coords  # lit at the vertices

        # soft sky light: brighter from above, never black below
        sky = 0.35 + 0.65 * (0.5 + 0.5 * (normals @ up))
        # hard key light from the upper left, in front: faces turned away are darker
        key_direction = _normalized((0.7 * up - 0.6 * right + 0.4 * back)[None, :])[0]
        key = np.clip(normals @ key_direction, 0.0, 1.0)
        # fill light from the lower right, opposite of the key
        fill_direction = _normalized((-0.4 * up + 0.6 * right + 0.3 * back)[None, :])[0]
        fill = np.clip(normals @ fill_direction, 0.0, 1.0)
        # spot light above the focus, in a cone that fades out between 12 and 50
        # degrees; wrapped (half) Lambert shading, soft over the cells
        spot_position = focus + 6.0 * a * up + 3.0 * a * back
        axis = _normalized((focus - spot_position)[None, :])[0]
        to_cell = _normalized(centers - spot_position[None, :])
        cone = _smoothstep(
            np.cos(np.radians(50.0)), np.cos(np.radians(12.0)), to_cell @ axis
        )
        spot = cone * (0.5 + 0.5 * np.clip(-_row_dots(normals, to_cell), 0.0, 1.0))
        # depth: elements behind the point the camera looks at fade
        distance = _row_norms(centers - camera[None, :])
        reference = np.linalg.norm(look_at - camera)
        fog = np.clip(1.0 - 0.5 * (distance - reference) / reference, 0.4, 1.0)

        lit = np.clip(
            self.SKY * sky + self.KEY * key + self.FILL * fill + self.SPOT * spot,
            0.0,
            1.0,
        )
        lit = np.where(has_normal, lit, 0.7) * fog
        self._values = self._kinds + np.where(self._lined, 0.5, 0.0) + self.SPAN * lit

    @property
    def has_envelope(self):
        return self._rx is not None

    def points(self):
        """All vertices of the elements, and the orbit, (num, 3)"""
        import numpy as np

        return np.concatenate((self._coords, self._p))

    def add_wires(self, radius, ring_spacing, camera=None):
        """The orbit and the envelope cage as thin opaque wires

        Ascent draws lines with a radius of 0.001 times their extent, too thick
        for close-ups: wires are thin boxes of any radius.

        Parameters
        ----------
        radius : float
            Half thickness of the wires [m].
        ring_spacing : float
            Least distance in s between the rings of the envelope cage [m].
        camera : tuple, optional
            Position of a camera close to the wires, e.g., inside the lattice [m], and
            the distance from it at which the wires are radius thick [m]. Their
            thickness then grows in proportion to the distance from the camera: they
            look as thin close by as far away.
        """
        import numpy as np

        if camera is None:

            def half_thickness(points):
                return np.full(len(points), radius)

        else:
            eye = np.asarray(camera[0], dtype=np.float64)
            per_distance = radius / camera[1]

            def half_thickness(points):
                return per_distance * _row_norms(points - eye[None, :])

        p, t, n, v, s = self._p, self._t, self._n, self._v, self._s
        thickness = half_thickness(p)
        parts = [(_box_sweep_varying(p, n, v, thickness, thickness), WIRE_ORBIT)]

        if self.has_envelope:
            rx, ry = self._rx, self._ry
            # longitudinal wires
            for k in range(_CAGE_WIRES):
                phi = 2.0 * np.pi * k / _CAGE_WIRES
                path = (
                    p
                    + (rx * np.cos(phi))[:, None] * n
                    + (ry * np.sin(phi))[:, None] * v
                )
                thickness = half_thickness(path)
                parts.append(
                    (
                        _box_sweep_varying(path, n, v, thickness, thickness),
                        WIRE_ENVELOPE,
                    )
                )
            # rings, at least ring_spacing apart from the start of the lattice and at
            # both of its ends
            rings = [0]
            for i in range(1, len(s)):
                if s[i] - s[rings[-1]] >= ring_spacing:
                    rings.append(i)
            if rings[-1] != len(s) - 1:
                rings.append(len(s) - 1)
            m = _RING_VERTICES
            phi = 2.0 * np.pi * np.arange(m + 1) / m  # closed
            cos, sin = np.cos(phi)[None, :, None], np.sin(phi)[None, :, None]
            # all rings in one sweep, (rings, m + 1, 3) as (rings * (m + 1), 3)
            rp, rn, rv = (
                p[rings][:, None, :],
                n[rings][:, None, :],
                v[rings][:, None, :],
            )
            around = (
                rp
                + rx[rings][:, None, None] * cos * rn
                + ry[rings][:, None, None] * sin * rv
            ).reshape(-1, 3)
            radial = (cos * rn + sin * rv).reshape(-1, 3)
            along = np.broadcast_to(t[rings][:, None, :], (len(rings), m + 1, 3))
            width = half_thickness(around)
            points, cells = _box_sweep_varying(
                around, along.reshape(-1, 3), radial, width, width
            )
            # without the cells that join the end of a ring to the start of the next
            joins = np.arange(len(cells)) % (m + 1) == m
            parts.append(((points, cells[~joins]), WIRE_ENVELOPE))

        coords, conn, values = [], [], []
        num_vertices = 0
        for (points, cells), value in parts:
            coords.append(points)
            conn.append(cells + num_vertices)
            values.append(np.full(len(cells), value))
            num_vertices += len(points)
        self._wires = (
            np.concatenate(coords),
            np.concatenate(conn),
            np.concatenate(values),
        )

    def add_anchor(self, position):
        """A point that the bounds of the surfaces contain, e.g., a camera in the lattice

        Viskores traces rays only for the pixels that the corners of the bounds of a
        surface cover, unless the camera is inside the bounds; corners behind the
        camera end up on the wrong side of the image, and near surfaces are cut off.
        A degenerate cell at the point, which no ray hits, extends the bounds of the
        surfaces of the opaque elements and of the wires.

        Parameters
        ----------
        position : array_like
            The point [m].
        """
        import numpy as np

        self._anchor = np.asarray(position, dtype=np.float64)

    def _anchored(self, points, cells, values, per_vertex):
        """Points, cells and values with a degenerate cell at the anchor"""
        import numpy as np

        if self._anchor is None:
            return points, cells, values
        corners = cells.shape[1]
        points = np.concatenate((points, np.repeat(self._anchor[None, :], corners, 0)))
        cells = np.concatenate(
            (cells, np.arange(corners)[None, :] + len(points) - corners)
        )
        # the value of the degenerate cell is never drawn
        extra = np.full(corners if per_vertex else 1, values[0])
        return points, cells, np.concatenate((values, extra))

    def blueprint(self):
        """The Blueprint domain, as a dict of numpy arrays and strings"""
        import numpy as np

        mesh = {"coordsets": {}, "topologies": {}, "fields": {}}
        if len(self._conn):
            coords, conn, values = self._coords, self._conn, self._values
            if self._surfaces:
                # surfaces: ray traced, see add_anchor
                coords, conn, values = self._anchored(coords, conn, values, True)
            _add_hexahedra(
                mesh,
                SURVEY_TOPOLOGY,
                coords,
                conn,
                shape="quad" if self._surfaces else "hex",
            )
            mesh["fields"][SURVEY_FIELD] = {
                "association": "vertex",
                "topology": SURVEY_TOPOLOGY,
                "values": np.asarray(values, dtype=np.float64),
            }
        if self._wires is not None:
            points, cells, values = self._anchored(*self._wires, False)
            _add_hexahedra(mesh, WIRES_TOPOLOGY, points, cells)
            mesh["fields"][WIRES_FIELD] = {
                "association": "element",
                "topology": WIRES_TOPOLOGY,
                "values": np.asarray(values, dtype=np.float64),
            }
        return mesh

    def diagonal(self):
        """Diagonal of the bounding box of the elements: the extent that the volume
        renderer refers its opacity to"""
        import numpy as np

        if not len(self._coords):
            return 0.0
        return float(
            np.linalg.norm(self._coords.max(axis=0) - self._coords.min(axis=0))
        )

    def samples(self, least=2000, most=4000):
        """Samples of the volume renderer: a few across the thinnest wall"""
        if not self._walls:
            return least
        thinnest = min(self._walls.values())
        if self._liner is not None:
            thinnest = min(thinnest, self._liner)
        return int(min(max(math.ceil(4.0 * self.diagonal() / thinnest), least), most))

    def color_table(self, opacity):
        """An Ascent color table of the elements: one control point per color index

        The volume renderer of VTK-h scales the opacity of the color table to the
        opacity across a tenth of the diagonal of the data. The alpha values here
        compensate for the wall thickness, so that each wall has about the given
        opacity, as far as the float32 alpha of the color table allows.
        """
        diagonal = self.diagonal()

        def alpha(target, thickness):
            if target <= 0.0 or thickness <= 0.0:
                return 0.0
            exponent = diagonal / (10.0 * thickness)
            # Viskores stores alpha as float32: 1 - alpha must stay a positive float32,
            # else it rounds to an opaque 1
            return min(1.0 - (1.0 - target) ** exponent, 1.0 - 2.0**-24)

        last = self.max_value
        control_points = []

        def point(value, rgb, a):
            control_points.append(
                {"type": "rgb", "position": value / last, "color": list(rgb)}
            )
            control_points.append(
                {"type": "alpha", "position": value / last, "alpha": a}
            )

        shade = _shade

        # without lighting, all surfaces are lit
        unlit = self.AMBIENT if self._lit else 1.0
        dark_tint = self.SHADOW_TINT if self._cinematic else (0.0, 0.0, 0.0)
        bright_tint = self.HIGHLIGHT_TINT if self._cinematic else (0.0, 0.0, 0.0)
        bright = 0.85 if self._cinematic else 1.0
        if self._cinematic:
            unlit = 0.25
        liner_shade = self.LINER_SHADE_LIT if self._lit else self.LINER_SHADE
        # cinematic: lit surfaces are denser than the ones in the shade, as if they
        # scattered the light
        dim, glow = (0.6, 1.8) if self._cinematic else (1.0, 1.0)
        point(0.0, (1.0, 1.0, 1.0), 0.0)
        for k, rgb in enumerate(self.colors, start=1):
            wall = self._walls[k]
            point(k, shade(rgb, unlit, dark_tint), alpha(dim * opacity, wall))
            point(
                k + self.SPAN,
                shade(rgb, bright, bright_tint),
                alpha(min(glow * opacity, 0.85), wall),
            )
            if self._liner is not None:
                # the lining of the bore: darker and denser, as a shaded inner wall
                liner = self._liner
                point(
                    k + 0.5,
                    shade(rgb, liner_shade * unlit, dark_tint),
                    alpha(dim * self.LINER_OPACITY, liner),
                )
                point(
                    k + 0.5 + self.SPAN,
                    shade(rgb, liner_shade * bright, bright_tint),
                    alpha(min(glow * self.LINER_OPACITY, 0.85), liner),
                )
        return {"control_points": control_points}

    def surface_value_range(self):
        """The range of the field values of opaque elements"""
        return 1.0, len(self.colors) + self.SPAN

    def surface_color_table(self):
        """An Ascent color table of opaque elements: per color, from dark to lit

        The field values hold the shading, see light and light_cinematic: the
        renderer does not shade the elements, it would show the flat faces.
        """
        lowest, highest = self.surface_value_range()
        control_points = []
        for k, rgb in enumerate(self.colors, start=1):
            if self._cinematic:
                dark = _shade(rgb, 0.25, self.SHADOW_TINT)
                bright = _shade(rgb, 1.0, self.HIGHLIGHT_TINT)
            elif self._lit:
                dark, bright = _shade(rgb, self.OPAQUE_AMBIENT), rgb
            else:
                dark = bright = rgb
            for value, color in ((k, dark), (k + self.SPAN, bright)):
                control_points.append(
                    {
                        "type": "rgb",
                        "position": (value - lowest) / (highest - lowest),
                        "color": list(color),
                    }
                )
        return {"control_points": control_points}

    def wires_color_table(self, orbit=(0.0, 0.0, 0.0)):
        """An Ascent color table of the wires: the orbit color, the envelope color"""
        return {
            "control_points": [
                {"type": "rgb", "position": WIRE_ORBIT, "color": list(orbit)},
                {
                    "type": "rgb",
                    "position": WIRE_ENVELOPE,
                    "color": list(ENVELOPE_COLOR),
                },
            ]
        }
