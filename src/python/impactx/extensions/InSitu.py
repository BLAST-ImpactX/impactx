"""
This file is part of ImpactX

Copyright 2026 ImpactX contributors
Authors: Axel Huebl, Chad Mitchell
License: BSD-3-Clause-LBNL
"""

import os
import shutil
import tempfile
import weakref

from ..impactx_pybind import CoordSystem

#: Visualization presets, see :py:func:`ascent_actions`
PRESETS = ("beam3d", "phasespace", "survey3d")

# plot units, as in ImpactXParticleContainer.plot_phasespace
_M2MM = 1.0e3
_RAD2MRAD = 1.0e3


def _my_rank():
    """MPI rank of this process in the communicator of ImpactX"""
    from ..impactx_pybind import amr

    return amr.ParallelDescriptor.MyProc()


def _host_accessible(pc):
    """Is the particle data of the beam accessible on the host (CPU)?"""
    from ..impactx_pybind import amr

    arena = pc.arena
    if arena is None:  # the default arena
        arena = amr.The_Arena()
    return arena.is_host_accessible


def _bin_range(vmin, vmax):
    """A non-empty bin range, also for a beam without spread in a coordinate"""
    vmin, vmax = float(vmin), float(vmax)
    if vmax > vmin:
        return vmin, vmax
    half = 0.5 * abs(vmin) if vmin != 0.0 else 1.0e-12
    return vmin - half, vmax + half


class _Visualization:
    """A visualization: the Ascent actions of a frame and what it shows"""

    def actions(self, pc=None, image_dir=None, element=None):
        """Ascent actions for the beam. MPI-collective.

        ``element`` is the element of the lattice that the beam is in or has just
        passed, if known.
        """
        raise NotImplementedError

    def geometry(self):
        """A Blueprint domain to publish next to the beam, or None"""
        return None

    def frame(self, info, renderer, frame_name):
        """The frame: what Ascent rendered"""
        return {"images": renderer.collect_images(info, frame_name)}


class _Beam3D(_Visualization):
    """Ascent renders the beam particles in x, y and t, colored by a particle component

    The extent of a beam in t and in x, y can differ by orders of magnitude, so
    by default, each axis is scaled to the extent of the beam on the first frame.
    """

    #: an oblique view on the beam
    DEFAULT_CAMERA = {"azimuth": 35.0, "elevation": 25.0}

    def __init__(
        self,
        color="momentum_t",
        image_size=(1024, 768),
        color_table=None,
        radius=None,
        camera=None,
        scale="auto",
    ):
        if not (scale is None or scale == "auto" or len(scale) == 3):
            raise ValueError(
                "beam3d: scale must be 'auto', None or three factors for x, y, t"
            )
        self.color = color
        self.image_size = image_size
        self.color_table = color_table
        self.radius = radius
        self.camera = camera
        self.scale = scale
        self._factors = None if scale == "auto" else scale

    def _scale_factors(self, pc):
        """Per axis, 1 / extent of the beam on the first frame. MPI-collective."""
        if self._factors is None:
            if pc is None:
                raise ValueError("beam3d: scale='auto' needs the beam.")
            rbc = pc.beam_moments()
            factors = []
            for name in ("x", "y", "t"):
                extent = float(rbc[f"max_{name}"]) - float(rbc[f"min_{name}"])
                factors.append(1.0 / extent if extent > 0.0 else 1.0)
            self._factors = factors
        return self._factors

    def actions(self, pc=None, image_dir=None, element=None):
        width, height = self.image_size
        plot = {"type": "pseudocolor", "field": self.color}
        if self.color_table is not None:
            plot["color_table"] = {"name": self.color_table}
        if self.radius is not None:
            plot["points"] = {"radius": self.radius}

        name = "beam3d" if image_dir is None else os.path.join(image_dir, "beam3d")
        render = {
            "image_name": name,
            "image_width": width,
            "image_height": height,
            "camera": self.DEFAULT_CAMERA if self.camera is None else self.camera,
        }

        actions = []
        if self.scale is not None:
            sx, sy, st = self._scale_factors(pc)
            actions.append(
                {
                    "action": "add_pipelines",
                    "pipelines": {
                        "beam3d": {
                            "f1": {
                                "type": "transform",
                                "params": {"scale": {"x": sx, "y": sy, "z": st}},
                            }
                        }
                    },
                }
            )
            plot["pipeline"] = "beam3d"
            # the axis values would be the scaled ones
            render["world_annotations"] = "false"

        actions.append(
            {
                "action": "add_scenes",
                "scenes": {
                    "beam3d": {"plots": {"p1": plot}, "renders": {"r1": render}}
                },
            }
        )
        return actions


class _OwnActions(_Visualization):
    """Own Ascent actions: the frames are the images that they render"""

    def __init__(self, actions):
        self._actions = actions

    def actions(self, pc=None, image_dir=None, element=None):
        return self._actions


class _PhaseSpace(_Visualization):
    """Ascent bins the x-px, y-py and t-pt projections, matplotlib draws them

    The particle charge is summed per bin over all MPI ranks, as shown by
    ImpactXParticleContainer.plot_phasespace. Ascent bins on the host (CPU):
    on GPUs, this needs the particles in managed memory.
    """

    def __init__(self, num_bins=50):
        if num_bins < 1:
            raise ValueError("phasespace: num_bins must be >= 1")
        self.num_bins = num_bins
        self._ranges = None
        self._rbc = None
        self._charge = None

    def actions(self, pc=None, image_dir=None, element=None):
        if pc is None:
            raise ValueError("ascent_actions: the 'phasespace' preset needs the beam.")
        if not _host_accessible(pc):
            raise RuntimeError(
                "phasespace: Ascent bins on the host (CPU), but the beam is in GPU "
                "device memory. Use managed memory (amrex.the_arena_is_managed = 1) "
                "or plot_phasespace()."
            )

        # bin ranges: MPI-collective
        rbc = pc.beam_moments()
        if rbc["charge_C"] == 0.0:
            raise RuntimeError(
                "phasespace: the beam has no charge. The in situ phase space shows the "
                "charge per bin; use plot_phasespace() to count macroparticles."
            )
        self._rbc = rbc
        self._charge = abs(pc.ref.charge)

        n = self.num_bins
        self._ranges = {}
        queries = {}
        for name in ("x", "y", "t"):
            rmin, rmax = _bin_range(rbc[f"min_{name}"], rbc[f"max_{name}"])
            pmin, pmax = _bin_range(rbc[f"min_p{name}"], rbc[f"max_p{name}"])
            self._ranges[name] = (rmin, rmax, pmin, pmax)
            # clamp: the particles on the upper end of a range go in the last bin
            expression = (
                "binning('weighting', 'sum', ["
                f"axis('position_{name}', num_bins={n}, "
                f"min_val={rmin:.17e}, max_val={rmax:.17e}, clamp=True), "
                f"axis('momentum_{name}', num_bins={n}, "
                f"min_val={pmin:.17e}, max_val={pmax:.17e}, clamp=True)])"
            )
            queries[f"phasespace_{name}"] = {
                "params": {"expression": expression, "name": f"phasespace_{name}"}
            }

        return [{"action": "add_queries", "queries": queries}]

    def frame(self, info, renderer, frame_name):
        if _my_rank() != 0:
            return {}

        import numpy as np

        n = self.num_bins
        hists = []
        for name in ("x", "y", "t"):
            # results by query name, then by cycle: the last cycle is this execute
            results = info["expressions"][f"phasespace_{name}"]
            result = list(results.values())[-1]
            values = np.asarray(result["attrs"]["value"]["value"], dtype=np.float64)
            # the first axis (position) varies fastest: rows along the position
            hist = values.reshape(n, n).T * self._charge

            rmin, rmax, pmin, pmax = self._ranges[name]
            p_unit = 1.0 if name == "t" else _RAD2MRAD
            r_edges = np.linspace(rmin * _M2MM, rmax * _M2MM, n + 1)
            p_edges = np.linspace(pmin * p_unit, pmax * p_unit, n + 1)
            hists.append((hist, r_edges, p_edges))

        return {"phasespace": (self._rbc, r"$|Q|$ [C/bin]", *hists)}


#: rms moments and plot labels of the particle components, for the color of the beam
_COMPONENT_MOMENTS = {
    "position_x": ("x", r"$\Delta x$ [m]"),
    "position_y": ("y", r"$\Delta y$ [m]"),
    "position_t": ("t", r"$\Delta ct$ [m]"),
    "momentum_x": ("px", r"$\Delta p_x$ [rad]"),
    "momentum_y": ("py", r"$\Delta p_y$ [rad]"),
    "momentum_t": ("pt", r"$\Delta p_t$ [$p_0\cdot c$]"),
}


def _draw_survey3d(ax, path, meta, legend=True, legend_ncols=5):
    """Draw a rendered survey into matplotlib axes, with legend and color bar"""
    import matplotlib.image as mpimg
    from matplotlib.patches import Patch

    ax.set_axis_off()
    ax.imshow(mpimg.imread(path))
    if legend and meta["legend"]:
        handles = [Patch(color=rgb, label=label) for label, rgb in meta["legend"]]
        ax.legend(
            handles=handles,
            bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
            loc="lower left",
            ncols=legend_ncols,
            mode="expand",
            borderaxespad=0.0,
        )
    if meta["colorbar"] is not None:
        import matplotlib as mpl

        label, vmin, vmax = meta["colorbar"]
        cax = ax.inset_axes([0.72, 0.12, 0.25, 0.03])
        cbar = mpl.colorbar.ColorbarBase(
            cax,
            cmap=mpl.colormaps["coolwarm"],
            norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax),
            orientation="horizontal",
            label=label,
            ticks=mpl.ticker.MaxNLocator(nbins=3),
            format=mpl.ticker.ScalarFormatter(useMathText=True),
        )
        cbar.formatter.set_powerlimits((-2, 2))
        cbar.update_ticks()
        if meta.get("dark"):
            # on the dark image
            cax.tick_params(colors="0.9")
            cax.xaxis.label.set_color("0.9")
            cax.xaxis.get_offset_text().set_color("0.9")
    if meta["note"] is not None:
        ax.text(
            0.01,
            0.02,
            meta["note"],
            transform=ax.transAxes,
            fontsize="small",
            color="0.9" if meta.get("dark") else "black",
        )


def _preset(preset, actions, kwargs):
    """The visualization of a preset or of own actions"""
    if actions is not None:
        if kwargs:
            raise TypeError(
                f"Options {sorted(kwargs)} are for presets, not for own actions."
            )
        return _OwnActions(actions)
    if preset == "beam3d":
        return _Beam3D(**kwargs)
    if preset == "phasespace":
        return _PhaseSpace(**kwargs)
    if preset == "survey3d":
        from .InSituSurvey3D import _Survey3D

        return _Survey3D(**kwargs)
    raise ValueError(f"Unknown preset '{preset}', choose one of {PRESETS}.")


def ascent_actions(preset, pc=None, **kwargs):
    """
    Ascent actions of a visualization preset.

    Use the result as a starting point for own Ascent actions.

    Parameters
    ----------
    preset : str
        ``"beam3d"``: Ascent renders the particles in x, y and t, colored by a
        particle component. Options: ``color`` (default: ``"momentum_t"``),
        ``image_size`` (pixels, default: ``(1024, 768)``), ``color_table`` (name),
        ``radius`` (of a point, in the scaled coordinates), ``camera`` (Ascent camera
        parameters, default: an oblique view) and ``scale``: ``"auto"`` (default)
        scales each axis by the inverse extent of the beam on the first frame, and
        renders no axes; ``None`` renders the beam in meters, with axes; or three
        factors for x, y and t.

        ``"phasespace"``: Ascent bins the x-px, y-py and t-pt phase space
        projections of the beam charge. Options: ``num_bins`` (per axis, default: 50).
        Ascent bins on the host (CPU): on GPUs, this needs the particles in
        managed memory.

        ``"survey3d"``: Ascent renders the lattice in the global frame of its survey,
        with the reference orbit, and with a beam, the beam particles and an
        envelope cage, see ``KnownElementsList.plot_survey(projection="3D")``. The
        actions refer to the geometry that ``Animation`` and ``plot_survey`` publish
        next to the beam. Options: ``lattice`` (required), ``ref`` (the reference
        particle, without a beam), ``envelope``, ``magnify``, ``element_size``,
        ``camera``, ``camera_smoothing``, ``image_size``, ``supersample``, ``palette``,
        ``opacity`` (of the element walls, default: 0.3; 1: opaque elements, closed
        around the bore, with drifts as beam pipes), ``color`` (particle
        component, default: ``"momentum_t"``), ``radius`` (of a particle [m]),
        ``samples`` (of the volume rendering), ``bore``, ``cutaway`` and
        ``lighting``.
    pc : ImpactXParticleContainer, optional
        The particle beam, needed for ``"phasespace"`` to set the bin ranges and for
        the beam in ``"survey3d"``. This is then MPI-collective.
    **kwargs
        Options of the preset.

    Returns
    -------
    list
        Ascent actions.

    Raises
    ------
    ValueError
        If the preset is unknown, or ``"phasespace"`` is requested without a beam.
    RuntimeError
        If ``"phasespace"`` is requested for a beam without charge, or for
        particles in GPU device memory.
    """
    return _preset(preset, None, kwargs).actions(pc)


class _Renderer:
    """An open Ascent instance, which collects the images it renders

    Ascent writes images on MPI rank 0. Each frame moves them to unique file
    names, so the Ascent actions can stay the same between frames and Ascent
    does not rebuild its pipelines.
    """

    def __init__(self, options=None, directory=None):
        from ..impactx_pybind import insitu

        self._rank = _my_rank()
        self.directory = None
        self.image_dir = None
        if self._rank == 0:
            if directory is None:
                self.directory = tempfile.mkdtemp(prefix="impactx-ascent-")
                weakref.finalize(
                    self, shutil.rmtree, self.directory, ignore_errors=True
                )
            else:
                self.directory = os.path.abspath(directory)
                os.makedirs(self.directory, exist_ok=True)
            # Ascent renders into this directory
            self.image_dir = os.path.join(self.directory, "ascent")
            os.makedirs(self.image_dir, exist_ok=True)

        self._ascent = insitu.Ascent(options if options is not None else {})

    def execute(self, pc, actions, step, period, geometry=None):
        """Publish the beam and a geometry, execute Ascent actions. MPI-collective."""
        if pc is not None and pc.coord_system != CoordSystem.s:
            raise RuntimeError(
                "Ascent: the particles must be at fixed s, "
                f"but are at {pc.coord_system}."
            )
        return self._ascent.execute(pc, actions, step, period, geometry)

    def collect_images(self, info, frame_name):
        """Move the images rendered by the last execute to unique names

        Returns the image files on MPI rank 0, else an empty list.
        """
        if self._rank != 0:
            return []

        images = info.get("images", [])
        if isinstance(images, dict):
            images = list(images.values())
        paths = []
        for k, image in enumerate(images):
            dst = os.path.join(self.directory, f"{frame_name}_{k:02d}.png")
            os.replace(image["image_name"], dst)
            paths.append(dst)
        return paths

    def close(self):
        """Close Ascent. MPI-collective. Rendered images stay."""
        self._ascent.close()


def _figure(frame, dpi=100):
    """A new, empty matplotlib figure that fits a frame"""
    import matplotlib.pyplot as plt

    if "phasespace" in frame:
        fig = plt.figure(figsize=(16, 4), constrained_layout=True)
    elif "survey3d" in frame:
        import matplotlib.image as mpimg

        height, width = mpimg.imread(frame["images"][0]).shape[:2]
        # room for the legend above the image
        fig = plt.figure(figsize=(width / dpi, 1.15 * height / dpi), dpi=dpi)
    else:
        import matplotlib.image as mpimg

        sizes = [mpimg.imread(p).shape for p in frame["images"]]
        width = sum(s[1] for s in sizes)
        height = max(s[0] for s in sizes)
        fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.canvas.manager.set_window_title("ImpactX: in situ")
    return fig


def _draw(fig, frame):
    """Draw a frame into an empty figure"""
    if "phasespace" in frame:
        from .ImpactXParticleContainer import _draw_phasespace

        _draw_phasespace(fig, *frame["phasespace"])
        return
    if "survey3d" in frame:
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0 / 1.15])
        _draw_survey3d(ax, frame["images"][0], frame["survey3d"])
        return

    import matplotlib.image as mpimg

    paths = frame["images"]
    axes = fig.subplots(1, len(paths), squeeze=False)[0]
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.0)
    for ax, path in zip(axes, paths):
        ax.set_axis_off()
        ax.imshow(mpimg.imread(path))


def _frame_rgb(frame, dpi):
    """A frame drawn into an RGB image, as a numpy array (height, width, 3)"""
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    fig = _figure(frame, dpi=dpi)
    try:
        _draw(fig, frame)
        # off screen, whichever backend is in use
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        return np.asarray(canvas.buffer_rgba())[..., :3].copy()
    finally:
        plt.close(fig)


def _has_data(frame):
    return bool(frame.get("images")) or "phasespace" in frame


def ix_pc_plot_ascent(
    self,
    preset="beam3d",
    *,
    actions=None,
    options=None,
    step=0,
    period=0,
    **kwargs,
):
    """
    Visualize the particle beam in situ with Ascent, as a matplotlib figure.

    ImpactX passes the particle arrays to Ascent without copies, relative to the
    reference particle at fixed s, over all MPI ranks.
    The coordinates are ``position_x``, ``position_y`` and ``position_t``, and every
    particle component is a field, e.g., ``momentum_t`` or ``weighting``.
    See :ref:`dataanalysis-insitu` for the copies that Ascent makes on GPUs.

    This is MPI-collective.

    Parameters
    ----------
    self : ImpactXParticleContainer
        The particle beam, at fixed s.
    preset : str, default="beam3d"
        A visualization preset, see :py:func:`impactx.insitu.ascent_actions`.
    actions : list, optional
        Ascent actions to execute instead of a preset. The figure shows the
        images they render.
    options : dict, optional
        Ascent options, see :py:class:`impactx.insitu.Ascent`.
    step : int, default=0
        Global step, the Ascent cycle.
    period : int, default=0
        Current period (turn or cycle).
    **kwargs
        Options of the preset.

    Returns
    -------
    A matplotlib figure.
    For MPI-parallel runs, the figure is only created on MPI rank 0.

    Raises
    ------
    RuntimeError
        If ImpactX is built without Ascent, Ascent fails, the particles are
        not at fixed s, or ``"phasespace"`` is requested for particles in GPU
        device memory.
    """
    vis = _preset(preset, actions, kwargs)

    renderer = _Renderer(options)
    try:
        actions = vis.actions(self, renderer.image_dir)
        info = renderer.execute(self, actions, step, period, vis.geometry())
        frame = vis.frame(info, renderer, "image")
    finally:
        renderer.close()

    if not _has_data(frame):
        return None
    fig = _figure(frame)
    _draw(fig, frame)
    return fig


class Animation:
    """
    Visualize the particle beam in situ with Ascent while it is tracked, and animate it.

    Each call renders one frame. Attach it to a tracking hook of the simulation,
    e.g., ``sim.hook["after_element"] = anim``, or call it with the simulation.

    The frames are kept on MPI rank 0. Show them in a window with :py:meth:`show`,
    save them as a GIF, an MP4 video, an HTML page or a directory of PNG images
    with :py:meth:`save`. In Jupyter, the animation displays itself.

    Ascent is opened on the first frame and stays open until :py:meth:`close`
    or until the simulation is finalized.

    Parameters
    ----------
    preset : str, default="beam3d"
        A visualization preset, see :py:func:`impactx.insitu.ascent_actions`.
    actions : list, optional
        Ascent actions to execute instead of a preset. The frames are the
        images they render.
    options : dict, optional
        Ascent options, see :py:class:`impactx.insitu.Ascent`.
    every : int, default=1
        Render every N-th call.
    directory : str, optional
        Where to keep rendered images. Default: a temporary directory that is
        removed with this object.
    **kwargs
        Options of the preset.
    """

    def __init__(
        self,
        preset="beam3d",
        *,
        actions=None,
        options=None,
        every=1,
        directory=None,
        **kwargs,
    ):
        if every < 1:
            raise ValueError("Animation: every must be >= 1")
        self._vis = _preset(preset, actions, kwargs)
        self._options = options
        self._every = every
        self._directory = directory
        self._calls = 0
        self._renderer = None
        #: Visualized frames, each with the step, the period, the position s and
        #: the time t (as c t) of the reference particle in m, and its data on MPI
        #: rank 0
        self.frames = []

    def __call__(self, sim):
        """Visualize the beam of a simulation, e.g., from a tracking hook.

        Renders every N-th call, see ``every``. This is MPI-collective.
        """
        self._calls += 1
        if (self._calls - 1) % self._every != 0:
            return
        self.render(
            sim.beam,
            step=sim.tracking_step,
            period=sim.tracking_period,
            element=sim.tracking_element,
        )

    def render(self, pc, step=0, period=0, element=None):
        """Visualize a particle beam as the next frame. This is MPI-collective.

        ``element`` is the lattice element that the beam is in or has just passed,
        for the ``"survey3d"`` preset, see :py:meth:`KnownElementsList.survey`.
        """
        if self._renderer is None:
            self._renderer = _Renderer(self._options, self._directory)
        actions = self._vis.actions(pc, self._renderer.image_dir, element=element)
        info = self._renderer.execute(pc, actions, step, period, self._vis.geometry())
        frame = self._vis.frame(info, self._renderer, f"frame_{len(self.frames):06d}")
        frame.update(step=step, period=period, s=pc.ref.s, t=pc.ref.t)
        self.frames.append(frame)

    def close(self):
        """Close Ascent. The frames stay. This is MPI-collective."""
        if self._renderer is not None:
            self._renderer.close()

    def _figure_animation(self, interval, dpi):
        """A matplotlib figure and the animation of the frames in it"""
        frames = [f for f in self.frames if _has_data(f)]
        if not frames:
            return None, None

        from matplotlib.animation import FuncAnimation

        fig = _figure(frames[0], dpi=dpi)

        def update(i):
            fig.clf()
            _draw(fig, frames[i])
            return []

        anim = FuncAnimation(
            fig, update, frames=len(frames), interval=interval, blit=False
        )
        return fig, anim

    def animation(self, interval=100, dpi=100):
        """
        The frames as a matplotlib animation.

        Parameters
        ----------
        interval : float, default=100
            Delay between frames in milliseconds.
        dpi : float, default=100
            Resolution of the figure: rendered images are shown at their pixel size.

        Returns
        -------
        matplotlib.animation.FuncAnimation, or None on MPI ranks other than 0
        and without frames.
        """
        _, anim = self._figure_animation(interval, dpi)
        return anim

    def show(self, interval=100):
        """Show the animation in a window (matplotlib). Only on MPI rank 0."""
        import matplotlib.pyplot as plt

        _, anim = self._figure_animation(interval, dpi=100)
        if anim is None:
            return
        plt.show()

    def save(self, filename, fps=10, dpi=100, pace="frame"):
        """
        Save the animation. Only on MPI rank 0.

        Parameters
        ----------
        filename : str
            Ending ``.gif`` (needs Pillow; 256 colors, shared by all frames and
            dithered), ``.webp`` (needs Pillow with WebP; full color), ``.mp4`` (needs
            FFmpeg; H.264, visually lossless) or ``.html``. Otherwise, a directory that
            receives the frames as PNG images.
        fps : float, default=10
            Frames per second, on average with ``pace="s"`` or ``"t"``.
        dpi : float, default=100
            Resolution of the figure: rendered images are shown at their pixel size.
        pace : str, default="frame"
            ``"frame"``: each frame is shown for the same time. ``"s"`` or ``"t"``:
            each frame is shown for a time that is proportional to the distance in s,
            or in the time t, of the reference particle to the next frame, so that the
            beam moves at a steady pace even if the frames are not evenly spaced.
            Needs a ``.gif`` or ``.webp`` file.
        """
        if pace not in ("frame", "s", "t"):
            raise ValueError("Animation: pace must be 'frame', 's' or 't'.")
        ext = os.path.splitext(filename)[1].lower()
        if pace != "frame" and ext not in (".gif", ".webp"):
            raise ValueError(
                "Animation: pace='s' or 't' needs a .gif or .webp file, which can show "
                "each frame for its own time."
            )
        if _my_rank() != 0:
            return

        import matplotlib.pyplot as plt

        if ext in (".gif", ".webp"):
            self._save_pillow(filename, ext, fps, dpi, pace)
            return
        if ext not in (".mp4", ".html"):
            os.makedirs(filename, exist_ok=True)
            for i, frame in enumerate(f for f in self.frames if _has_data(f)):
                fig = _figure(frame, dpi=dpi)
                try:
                    _draw(fig, frame)
                    fig.savefig(os.path.join(filename, f"frame_{i:06d}.png"), dpi=dpi)
                finally:
                    plt.close(fig)
            return

        fig, anim = self._figure_animation(1000.0 / fps, dpi)
        if anim is None:
            return
        try:
            if ext == ".html":
                with open(filename, "w") as f:
                    f.write(anim.to_jshtml(fps=fps))
            else:
                from matplotlib.animation import FFMpegWriter

                # H.264, visually lossless: thin wires and small particles stay sharp
                writer = FFMpegWriter(
                    fps=fps, codec="h264", extra_args=["-crf", "18", "-preset", "slow"]
                )
                anim.save(filename, writer=writer, dpi=dpi)
        finally:
            plt.close(fig)

    def _durations(self, frames, fps, pace):
        """Display time of each frame in ms, see save()"""
        import numpy as np

        mean = 1000.0 / fps
        if pace == "frame" or len(frames) < 2:
            return [int(round(mean))] * len(frames)
        position = np.array([frame[pace] for frame in frames], dtype=np.float64)
        step = np.abs(np.diff(position))
        if not np.any(step > 0.0):
            return [int(round(mean))] * len(frames)
        # the last frame as long as an average one
        step = np.append(step, step.mean())
        durations = mean * step / step.mean()
        # frames at the same position, e.g., around a zero-length element, flash briefly
        return [int(round(d)) for d in np.maximum(durations, 20.0)]

    def _save_pillow(self, filename, ext, fps, dpi, pace="frame"):
        """Save the frames as an animated GIF or WebP image, with Pillow

        A GIF has 256 colors: one palette for all frames, from all frames, and
        Floyd-Steinberg dithering keep smooth gradients, e.g., of dark backgrounds,
        free of bands and flicker.
        """
        from PIL import Image, features

        if ext == ".webp" and not features.check("webp"):
            raise RuntimeError(
                "Animation: saving .webp needs Pillow with WebP support."
            )
        frames = [frame for frame in self.frames if _has_data(frame)]
        images = [Image.fromarray(_frame_rgb(frame, dpi)) for frame in frames]
        if not images:
            return
        duration = self._durations(frames, fps, pace)
        if ext == ".webp":
            images[0].save(
                filename,
                save_all=True,
                append_images=images[1:],
                duration=duration,
                loop=0,
                quality=92,
                method=4,
            )
            return

        try:
            mediancut, floyd_steinberg = (
                Image.Quantize.MEDIANCUT,
                Image.Dither.FLOYDSTEINBERG,
            )
        except AttributeError:  # Pillow < 9.1
            mediancut, floyd_steinberg = Image.MEDIANCUT, Image.FLOYDSTEINBERG
        # one palette for all frames, from up to 8 of them at full resolution: thin
        # features, e.g., the color bar, need their pixels to get their colors
        samples = images[:: max(1, len(images) // 8)][:8]
        montage = Image.new("RGB", (samples[0].width, sum(im.height for im in samples)))
        y = 0
        for im in samples:
            montage.paste(im, (0, y))
            y += im.height
        palette = montage.quantize(colors=256, method=mediancut)
        frames = [im.quantize(palette=palette, dither=floyd_steinberg) for im in images]
        frames[0].save(
            filename,
            save_all=True,
            append_images=frames[1:],
            duration=duration,
            loop=0,
            optimize=False,
        )

    def _repr_html_(self):
        """Jupyter: an interactive animation of the frames"""
        import matplotlib.pyplot as plt

        fig, anim = self._figure_animation(interval=100, dpi=100)
        if anim is None:
            return None
        try:
            return anim.to_jshtml()
        finally:
            # only show the animation, not also the figure of its first frame
            plt.close(fig)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def register_InSitu_extension(insitu, ixpc):
    """In situ visualization helpers"""
    insitu.PRESETS = PRESETS
    insitu.ascent_actions = ascent_actions
    insitu.Animation = Animation
    Animation.__module__ = insitu.__name__

    ixpc.plot_ascent = ix_pc_plot_ascent
