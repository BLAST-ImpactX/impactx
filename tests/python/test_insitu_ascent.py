#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

"""In situ visualization of the beam with Ascent."""

import importlib.util
import os

import numpy as np
import pytest

from impactx import (
    Config,
    CoordSystem,
    ImpactX,
    amr,
    coordinate_transformation,
    distribution,
    elements,
    insitu,
)

needs_ascent = pytest.mark.skipif(
    not Config.have_ascent, reason="ImpactX is built without Ascent"
)


def on_root():
    """Figures, images and frames are only created on MPI rank 0"""
    return amr.ParallelDescriptor.MyProc() == 0


def skip_unless_beam_on_host(sim):
    """Ascent bins on the host (CPU): on GPUs, the beam must be in managed memory"""
    arena = sim.beam.arena or amr.The_Arena()
    if not arena.is_host_accessible:
        pytest.skip("Ascent bins on the host, the beam is in GPU device memory")


def make_sim(charge_C=1.0e-9):
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.verbose = 0
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(100.0)
    distr = distribution.Waterbag(
        lambdaX=1.0e-4,
        lambdaY=2.0e-4,
        lambdaT=1.0e-3,
        lambdaPx=1.0e-5,
        lambdaPy=2.0e-5,
        lambdaPt=1.0e-3,
    )
    sim.add_particles(charge_C, distr, 2000)
    return sim


@pytest.fixture
def sim():
    sim = make_sim()
    yield sim
    sim.finalize()


def panels(fig):
    """The axes of a figure that show an image"""
    return [ax for ax in fig.axes if ax.images]


def image(panel):
    return np.asarray(panel.images[0].get_array())


@pytest.mark.skipif(Config.have_ascent, reason="ImpactX is built with Ascent")
def test_without_ascent():
    with pytest.raises(RuntimeError, match="without Ascent"):
        insitu.Ascent()


@needs_ascent
def test_phasespace_matches_plot_phasespace(sim):
    """Ascent bins the charge per bin like plot_phasespace does"""
    skip_unless_beam_on_host(sim)
    fig_ref = sim.beam.plot_phasespace(num_bins=20)
    fig = sim.beam.plot_ascent("phasespace", num_bins=20)
    if not on_root():
        assert fig is None
        return

    for p_ref, p in zip(panels(fig_ref), panels(fig), strict=True):
        h_ref = p_ref.images[0].get_array().filled(0.0)
        h = p.images[0].get_array().filled(0.0)
        np.testing.assert_allclose(h, h_ref, rtol=1.0e-12, atol=1.0e-12 * h_ref.max())
        np.testing.assert_allclose(
            p.images[0].get_extent(), p_ref.images[0].get_extent(), rtol=1.0e-12
        )
        assert p.images[0].colorbar.ax.get_ylabel() == "$|Q|$ [C/bin]"
        # all of the beam charge is binned
        assert h.sum() == pytest.approx(1.0e-9, rel=1.0e-12)


@needs_ascent
def test_phasespace_without_charge():
    """The in situ phase space shows charge, a beam without charge raises"""
    sim = make_sim(charge_C=0.0)
    try:
        skip_unless_beam_on_host(sim)
        with pytest.raises(RuntimeError, match="no charge"):
            sim.beam.plot_ascent("phasespace")
    finally:
        sim.finalize()


@needs_ascent
def test_phasespace_needs_host_memory():
    """Ascent bins on the host: a beam in GPU device memory raises"""

    class DeviceMemory:
        is_host_accessible = False

    class BeamInDeviceMemory:
        arena = DeviceMemory()

    with pytest.raises(RuntimeError, match="device memory"):
        insitu.ascent_actions("phasespace", BeamInDeviceMemory())


@needs_ascent
@pytest.mark.parametrize("scale", ["auto", None, (1.0e3, 1.0e3, 1.0e2)])
def test_beam3d(sim, scale):
    fig = sim.beam.plot_ascent("beam3d", image_size=(160, 120), scale=scale)
    if not on_root():
        assert fig is None
        return

    (panel,) = panels(fig)
    pixels = image(panel)
    assert pixels.shape[:2] == (120, 160)
    # the particles are drawn on the background
    assert np.ptp(pixels[..., :3]) > 0.0


@needs_ascent
def test_own_actions(sim):
    """Own actions: the figure shows every image they render"""
    renders = {
        f"r{i}": {
            "image_name": f"view{i}",
            "image_width": 64,
            "image_height": 48,
            "camera": {"azimuth": 45.0 * i},
        }
        for i in range(2)
    }
    actions = [
        {
            "action": "add_scenes",
            "scenes": {
                "s1": {
                    "plots": {"p1": {"type": "pseudocolor", "field": "weighting"}},
                    "renders": renders,
                }
            },
        }
    ]
    fig = sim.beam.plot_ascent(actions=actions)
    if on_root():
        assert len(panels(fig)) == 2
        # the images are moved out of the current directory
        assert not os.path.exists("view0.png")
        assert not os.path.exists("view1.png")

    with pytest.raises(TypeError, match="presets"):
        sim.beam.plot_ascent(actions=actions, image_size=(64, 48))


@needs_ascent
def test_unknown_preset(sim):
    with pytest.raises(ValueError, match="Unknown preset"):
        sim.beam.plot_ascent("nope")


@needs_ascent
def test_actions_file_in_current_directory_is_ignored(sim):
    """Ascent would replace the actions with an ascent_actions.yaml in the cwd"""
    with open("ascent_actions.yaml", "w") as f:
        f.write('-\n  action: "this_action_does_not_exist"\n')
    fig = sim.beam.plot_ascent("beam3d", image_size=(64, 48))
    if on_root():
        assert len(panels(fig)) == 1


@needs_ascent
def test_errors_raise(sim):
    """Errors in Ascent reach Python, instead of being printed"""
    # clamp is a bool
    expression = "binning('weighting', 'sum', [axis('position_x', clamp=1)])"
    actions = [
        {
            "action": "add_queries",
            "queries": {"q1": {"params": {"expression": expression, "name": "broken"}}},
        }
    ]
    with pytest.raises(RuntimeError, match="axis"):
        sim.beam.plot_ascent(actions=actions)


def cube(as_lists=False):
    """A Blueprint domain of one hexahedron, with a vertex field"""
    corners = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=np.float64,
    )

    def array(values):
        return values.tolist() if as_lists else values

    return {
        "coordsets": {
            "cube_coords": {
                "type": "explicit",
                "values": {axis: array(corners[:, k]) for k, axis in enumerate("xyz")},
            }
        },
        "topologies": {
            "cube": {
                "type": "unstructured",
                "coordset": "cube_coords",
                "elements": {"shape": "hex", "connectivity": array(np.arange(8))},
            }
        },
        "fields": {
            "cube_value": {
                "association": "vertex",
                "topology": "cube",
                "values": array(np.arange(8, dtype=np.float32)),
            }
        },
    }


@needs_ascent
def test_execute_geometry(tmp_path):
    """A geometry as a dict of numpy arrays, of lists, or as JSON text"""
    import json

    with insitu.Ascent() as ascent:
        geometries = {
            "arrays": cube(),
            "lists": cube(as_lists=True),
            "json": json.dumps(cube(as_lists=True)),
        }
        for name, geometry in geometries.items():
            path = os.path.join(tmp_path, f"cube_{name}")
            actions = [
                {
                    "action": "add_scenes",
                    "scenes": {
                        "s1": {
                            "plots": {
                                "p1": {"type": "pseudocolor", "field": "cube_value"}
                            },
                            "renders": {
                                "r1": {
                                    "image_name": path,
                                    "image_width": 32,
                                    "image_height": 32,
                                }
                            },
                        }
                    },
                }
            ]
            ascent.execute(None, actions, geometry=geometry)
            if on_root():
                assert os.path.isfile(path + ".png")

        # Conduit has no booleans, unsigned or object arrays
        bad = cube()
        bad["fields"]["cube_value"]["values"] = np.arange(8, dtype=np.uint8)
        with pytest.raises(TypeError, match="uint8"):
            ascent.execute(None, [], geometry=bad)
        bad["fields"]["cube_value"]["values"] = [True] * 8
        with pytest.raises(TypeError, match="bool"):
            ascent.execute(None, [], geometry=bad)
        bad = cube()
        bad["fields"]["cube_value"]["flag"] = True
        with pytest.raises(TypeError, match="booleans"):
            ascent.execute(None, [], geometry=bad)
        bad = cube()
        bad["fields"]["cube_value"]["values"] = [object()] * 8
        with pytest.raises(TypeError, match="object"):
            ascent.execute(None, [], geometry=bad)
        bad = cube()
        bad["fields"][1] = bad["fields"].pop("cube_value")
        with pytest.raises(TypeError, match="keys"):
            ascent.execute(None, [], geometry=bad)
        # without a beam, the geometry is needed
        with pytest.raises(ValueError, match="needs a beam or a geometry"):
            ascent.execute(None, [], geometry={})


@needs_ascent
def test_fixed_t_raises(sim):
    coordinate_transformation(sim.beam, CoordSystem.t)
    try:
        with pytest.raises(RuntimeError, match="fixed s"):
            sim.beam.plot_ascent("beam3d")
    finally:
        coordinate_transformation(sim.beam, CoordSystem.s)


@needs_ascent
@pytest.mark.parametrize("preset", ["beam3d", "phasespace"])
def test_animation(sim, tmp_path, preset):
    if preset == "phasespace":
        skip_unless_beam_on_host(sim)
    options = {"image_size": (80, 60)} if preset == "beam3d" else {"num_bins": 8}
    anim = insitu.Animation(preset, every=2, **options)
    sim.hook["after_element"] = anim

    sim.lattice.extend([elements.Drift(name=f"d{i}", ds=0.5) for i in range(5)])
    sim.track_particles()

    # every 2nd of 5 elements
    assert len(anim.frames) == 3
    s = [frame["s"] for frame in anim.frames]
    assert s == pytest.approx([0.5, 1.5, 2.5])
    steps = [frame["step"] for frame in anim.frames]
    assert steps == sorted(steps)

    anim.save(str(tmp_path / "frames"))
    anim.save(str(tmp_path / "anim.html"))
    has_pil = importlib.util.find_spec("PIL") is not None
    has_webp = False
    if has_pil:
        from PIL import features

        has_webp = features.check("webp")
        anim.save(str(tmp_path / "anim.gif"), fps=4)
        if has_webp:
            anim.save(str(tmp_path / "anim.webp"), fps=4)

    if on_root():
        assert sorted(os.listdir(tmp_path / "frames")) == [
            f"frame_{i:06d}.png" for i in range(3)
        ]
        assert "<img" in (tmp_path / "anim.html").read_text()
        assert anim._repr_html_() is not None
        if has_pil:
            from PIL import Image

            with Image.open(tmp_path / "anim.gif") as gif:
                assert gif.n_frames == 3
        if has_webp:
            with Image.open(tmp_path / "anim.webp") as webp:
                assert webp.n_frames == 3
    else:
        assert anim.animation() is None


@needs_ascent
def test_animation_closed_by_finalize():
    """Ascent is closed with the simulation, the frames stay"""
    sim = make_sim()
    anim = insitu.Animation("beam3d", image_size=(64, 48))
    anim(sim)
    sim.finalize()

    assert len(anim.frames) == 1
    assert not anim._renderer._ascent.is_open
    if on_root():
        assert "<img" in anim._repr_html_()


def survey_lattice():
    """A quadrupole and a sector bend, with drifts"""
    return [
        elements.Quad(name="qf", ds=0.3, k=2.0, nslice=2),
        elements.Drift(name="d1", ds=0.5),
        elements.Sbend(name="b1", ds=1.0, rc=4.0, nslice=2),
        elements.Drift(name="d2", ds=0.5),
    ]


@needs_ascent
@pytest.mark.parametrize(
    "lighting, bore, cutaway, opacity, camera",
    [
        (None, "round", False, 0.3, "overview"),
        ("camera", "square", False, 0.3, "overview"),
        ("cinematic", "round", False, 0.3, "overview"),
        (None, None, True, 0.3, "overview"),
        ("cinematic", "round", True, 0.3, "overview"),
        ("camera", None, False, 1.0, "overview"),
        ("cinematic", "round", False, 0.3, "racecar"),
        ("cinematic", "round", True, 0.3, "racecar"),
        ("cinematic", "round", False, 1.0, "racecar"),
    ],
)
def test_plot_survey_3d(sim, lighting, bore, cutaway, opacity, camera):
    """Ascent renders the lattice, the orbit, the envelope and the beam"""
    import matplotlib.pyplot as plt

    sim.lattice.extend(survey_lattice())
    ax = sim.lattice.plot_survey(
        projection="3D",
        beam=sim.beam,
        image_size=(160, 90),
        lighting=lighting,
        bore=bore,
        cutaway=cutaway,
        opacity=opacity,
        camera=camera,
    )
    if on_root():
        rendered = image(ax)
        assert rendered.shape[:2] == (90, 160)
        # not only the background
        assert np.ptp(rendered[..., :3]) > 0.0
        plt.close(ax.figure)
    else:
        assert ax is None


@needs_ascent
def test_plot_survey_3d_lattice_only(sim):
    """Without a beam: the lattice and the orbit of a reference particle"""
    import matplotlib.pyplot as plt

    sim.lattice.extend(survey_lattice())
    ax = sim.lattice.plot_survey(
        sim.beam.ref.copy(), projection="3D", image_size=(160, 90)
    )
    if on_root():
        assert np.ptp(image(ax)[..., :3]) > 0.0
        plt.close(ax.figure)


@needs_ascent
def test_plot_survey_3d_errors(sim):
    sim.lattice.extend(survey_lattice())
    lattice = sim.lattice
    with pytest.raises(ValueError, match="lighting"):
        lattice.plot_survey(projection="3D", beam=sim.beam, lighting="sunny")
    with pytest.raises(ValueError, match="bore"):
        lattice.plot_survey(projection="3D", beam=sim.beam, bore="hexagonal")
    with pytest.raises(ValueError, match="camera must be"):
        lattice.plot_survey(projection="3D", beam=sim.beam, camera="helicopter")
    for smoothing in ("wide", -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="camera_smoothing"):
            lattice.plot_survey(
                projection="3D",
                beam=sim.beam,
                camera="beam",
                camera_smoothing=smoothing,
            )
    for opacity in ("dense", -0.1, 1.5, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="opacity"):
            lattice.plot_survey(projection="3D", beam=sim.beam, opacity=opacity)
    with pytest.raises(ValueError, match="projection"):
        lattice.plot_survey(projection="XY")
    with pytest.raises(ValueError, match="not both"):
        lattice.plot_survey(sim.beam.ref.copy(), projection="3D", beam=sim.beam)


@needs_ascent
@pytest.mark.parametrize("camera", ["beam", "racecar"])
def test_animation_survey3d(sim, camera):
    """The survey follows the beam through the lattice, element by element"""
    sim.lattice.extend(survey_lattice())
    anim = insitu.Animation(
        "survey3d",
        lattice=sim.lattice,
        camera=camera,
        lighting="cinematic",
        image_size=(80, 45),
    )
    sim.hook["after_element"] = anim
    sim.track_particles()

    assert len(anim.frames) == len(survey_lattice())
    if on_root():
        assert all(len(frame["images"]) == 1 for frame in anim.frames)
        assert anim._repr_html_() is not None


def test_animation_pace():
    """With pace="s" or "t", frames are shown for times proportional to their spacing"""
    anim = insitu.Animation("beam3d")
    frames = [{"s": 0.0, "t": 0.0}, {"s": 1.0, "t": 2.0}, {"s": 3.0, "t": 3.0}]
    # spacing 1, 2, and the mean 1.5 for the last frame; 10 frames per second on average
    assert anim._durations(frames, 10.0, "s") == [67, 133, 100]
    assert anim._durations(frames, 10.0, "t") == [133, 67, 100]
    assert anim._durations(frames, 10.0, "frame") == [100, 100, 100]
    with pytest.raises(ValueError, match="pace"):
        anim.save("anim.gif", pace="z")
    with pytest.raises(ValueError, match="gif or .webp"):
        anim.save("anim.mp4", pace="s")
