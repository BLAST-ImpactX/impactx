#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

from pathlib import Path

import pytest

from impactx import ImpactX, distribution, elements


def test_diag_file_prefix():
    """
    This tests that diagnostics can be written under a custom file prefix.
    """
    sim = ImpactX()

    assert sim.diag_file_prefix == "diags"
    sim.diag_file_prefix = "initial_diags"
    assert sim.diag_file_prefix == "initial_diags"

    sim.particle_shape = 2
    sim.slice_step_diagnostics = True
    sim.init_grids()

    sim.diag_file_prefix = "custom_diags/"
    assert sim.diag_file_prefix == "custom_diags"

    # init particle beam
    kin_energy_MeV = 2.0e3
    bunch_charge_C = 1.0e-9
    npart = 1000

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

    distr = distribution.Waterbag(
        lambdaX=3.9984884770e-5,
        lambdaY=3.9984884770e-5,
        lambdaT=1.0e-3,
        lambdaPx=2.6623538760e-5,
        lambdaPy=2.6623538760e-5,
        lambdaPt=2.0e-3,
        muxpx=-0.846574929020762,
        muypy=0.846574929020762,
        mutpt=0.0,
    )
    sim.add_particles(bunch_charge_C, distr, npart)

    sim.lattice.append(elements.Drift(name="d1", ds=0.25))
    sim.track_particles()

    sim.finalize()

    custom_diags = Path("custom_diags")
    assert list(custom_diags.glob("ref_particle.*"))
    assert list(custom_diags.glob("reduced_beam_characteristics.*"))
    assert list(custom_diags.glob("ref_particle_final.*"))
    assert list(custom_diags.glob("reduced_beam_characteristics_final.*"))
    assert not list(Path("initial_diags").glob("ref_particle.*"))
    assert not list(Path("initial_diags").glob("reduced_beam_characteristics.*"))
    assert not Path("diags").exists()


@pytest.mark.parametrize("file_prefix", ["", "."])
def test_diag_file_prefix_current_directory(file_prefix):
    """
    This tests that an empty or "." diagnostic prefix writes to the current directory.
    """
    sim = ImpactX()

    sim.diag_file_prefix = file_prefix
    assert sim.diag_file_prefix == file_prefix

    Path("keep.txt").write_text("keep")

    sim.particle_shape = 2
    sim.slice_step_diagnostics = True
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    sim.lattice.append(elements.Drift(name="d1", ds=0.25))
    sim.track_reference(ref)

    sim.finalize()

    assert Path("keep.txt").read_text() == "keep"
    assert list(Path(".").glob("ref_particle.*"))
    assert list(Path(".").glob("ref_particle_final.*"))
    assert not Path("diags").exists()


def test_diag_file_prefix_sibling_directory_cleaned():
    """
    This tests that a sibling diagnostic prefix can be moved aside.
    """
    sibling_diags = Path("..") / f"{Path.cwd().name}_sibling_diags"
    sibling_diags.mkdir()
    (sibling_diags / "stale.txt").write_text("stale")

    sim = ImpactX()

    sim.diag_file_prefix = f"../{sibling_diags.name}/"
    assert sim.diag_file_prefix == f"../{sibling_diags.name}"

    sim.particle_shape = 2
    sim.init_grids()

    sim.finalize()

    assert sibling_diags.exists()
    assert not (sibling_diags / "stale.txt").exists()
    assert list(Path("..").glob(f"{sibling_diags.name}.old.*"))


def test_diag_file_prefix_internal_parent_directory_cleaned():
    """
    This tests that an internal ".." component still allows cleaning after normalization.
    """
    old_file_prefix = Path("external") / "diags"
    old_file_prefix.mkdir(parents=True)
    (old_file_prefix / "stale.txt").write_text("stale")

    sim = ImpactX()

    sim.diag_file_prefix = "external/abc/../diags/"
    assert sim.diag_file_prefix == "external/diags"

    sim.particle_shape = 2
    sim.init_grids()

    sim.finalize()

    assert old_file_prefix.exists()
    assert not (old_file_prefix / "stale.txt").exists()
    assert list(Path("external").glob("diags.old.*"))


def test_diag_file_prefix_parent_path_into_current_directory_cleaned():
    """
    This tests that a parent-relative path resolving back into CWD can be cleaned.
    """
    old_file_prefix = Path("diags")
    old_file_prefix.mkdir()
    (old_file_prefix / "stale.txt").write_text("stale")

    sim = ImpactX()

    sim.diag_file_prefix = f"../{Path.cwd().name}/diags/"
    assert sim.diag_file_prefix == f"../{Path.cwd().name}/diags"

    sim.particle_shape = 2
    sim.init_grids()

    sim.finalize()

    assert old_file_prefix.exists()
    assert not (old_file_prefix / "stale.txt").exists()
    assert list(Path(".").glob("diags.old.*"))
