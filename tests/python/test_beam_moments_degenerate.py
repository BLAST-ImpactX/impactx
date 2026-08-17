#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math
import os
import subprocess
import sys

import pytest

from impactx import ImpactX, create_envelope, distribution

# beam moments that stay defined without a beam: the reference path length and
# the beam charge (of which there is none)
ALWAYS_DEFINED = ("s", "charge_C")


def waterbag(lambda_pt=2.0e-3):
    """A 2 GeV electron beam, optionally without any energy spread."""
    return distribution.Waterbag(
        lambdaX=3.9984884770e-5,
        lambdaY=3.9984884770e-5,
        lambdaT=1.0e-3,
        lambdaPx=2.6623538760e-5,
        lambdaPy=2.6623538760e-5,
        lambdaPt=lambda_pt,
        muxpx=-0.846574929020762,
        muypy=0.846574929020762,
        mutpt=0.0,
    )


def test_beam_moments_without_beam():
    """
    An empty particle container has no moments to report.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    assert sim.beam.total_number_of_particles() == 0
    moments = sim.beam.beam_moments()

    undefined = {k: v for k, v in moments.items() if k not in ALWAYS_DEFINED}
    assert undefined
    assert all(math.isnan(v) for v in undefined.values())
    assert moments["charge_C"] == 0.0

    sim.finalize()


def test_beam_moments_without_reference_particle():
    """
    A beam that a Source element only loads later on has neither particles nor a
    reference particle when the initial diagnostics run. beta*gamma is not a real
    number for such a reference particle, which must not leak into the moments.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    assert sim.beam.ref.kin_energy_MeV == 0.0
    moments = sim.beam.beam_moments()

    undefined = {k: v for k, v in moments.items() if k not in ALWAYS_DEFINED}
    assert undefined
    assert all(math.isnan(v) for v in undefined.values())
    assert moments["charge_C"] == 0.0

    sim.finalize()


def test_beam_moments_without_charge():
    """
    Massless "test" particles carry no weight, so the beam has no weighted moments.
    Their extrema are unweighted, though, and stay defined.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(0.0, waterbag(), 1000)

    moments = sim.beam.beam_moments()

    assert moments["charge_C"] == 0.0
    assert math.isnan(moments["mean_x"])
    assert math.isnan(moments["sigma_x"])
    assert math.isnan(moments["emittance_x"])

    for coordinate in ("x", "y", "t", "px", "py", "pt"):
        assert moments[f"min_{coordinate}"] < moments[f"max_{coordinate}"]
        # deprecated spellings of the same two entries
        assert moments[f"{coordinate}_min"] == moments[f"min_{coordinate}"]
        assert moments[f"{coordinate}_max"] == moments[f"max_{coordinate}"]

    sim.finalize()


def test_beam_moments_single_particle():
    """
    A single particle has vanishing rms emittances, for which the Courant-Snyder
    (Twiss) functions are undefined.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(1.0e-9, waterbag(), 1)

    moments = sim.beam.beam_moments()

    # the beam is a point in phase space: no extent, no emittance
    for key in ("sigma_x", "sigma_y", "sigma_t", "emittance_x", "emittance_y"):
        assert moments[key] == 0.0

    # ... so its Twiss functions do not exist, in any plane
    for key in ("alpha_x", "alpha_y", "alpha_t", "beta_x", "beta_y", "beta_t"):
        assert math.isnan(moments[key])

    # its position is still perfectly well defined
    assert math.isfinite(moments["mean_x"])
    assert moments["mean_x"] == moments["min_x"] == moments["max_x"]

    sim.finalize()


def test_beam_moments_without_energy_spread():
    """
    A beam without energy spread has no longitudinal emittance, and hence no
    longitudinal Twiss functions. The transverse planes stay unaffected.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    sim.add_particles(1.0e-9, waterbag(lambda_pt=0.0), 1000)

    moments = sim.beam.beam_moments()

    assert moments["sigma_pt"] == 0.0
    assert moments["emittance_t"] == 0.0
    assert math.isnan(moments["beta_t"])
    assert math.isnan(moments["alpha_t"])

    assert moments["sigma_t"] > 0.0
    assert moments["emittance_x"] > 0.0
    assert math.isfinite(moments["beta_x"])
    assert math.isfinite(moments["alpha_x"])

    sim.finalize()


def test_envelope_beam_moments_without_energy_spread():
    """
    The covariance-matrix (envelope) moments treat a vanishing emittance the same
    way as the particle moments do.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(2.0e3)

    envelope = create_envelope(waterbag(lambda_pt=0.0), 1.0e-9)
    moments = envelope.beam_moments(ref)

    assert moments["emittance_t"] == 0.0
    assert math.isnan(moments["beta_t"])
    assert math.isnan(moments["alpha_t"])

    assert moments["emittance_x"] > 0.0
    assert math.isfinite(moments["beta_x"])
    assert math.isfinite(moments["alpha_x"])

    sim.finalize()


# A stand-alone run of the degenerate cases above, with AMReX trapping the
# floating-point exceptions that they used to raise. It has to run in its own
# process: a trap raises SIGFPE, which no test could survive in-process.
_FPE_TRAP_PROGRAM = """
from mpi4py import MPI  # noqa: F401  (keep MPI alive across ImpactX instances)

import amrex.space3d as amr

amr.initialize(
    [
        "amrex.verbose=0",
        "tiny_profiler.enabled=0",
        "amrex.signal_handling=1",
        "amrex.handle_sigfpe=1",
        "amrex.fpe_trap_invalid=1",
        "amrex.fpe_trap_zero=1",
    ]
)

from impactx import ImpactX, distribution


def beam(lambda_pt):
    return distribution.Waterbag(
        lambdaX=3.9984884770e-5,
        lambdaY=3.9984884770e-5,
        lambdaT=1.0e-3,
        lambdaPx=2.6623538760e-5,
        lambdaPy=2.6623538760e-5,
        lambdaPt=lambda_pt,
    )


def moments(npart, bunch_charge_C=1.0e-9, lambda_pt=2.0e-3, reference_particle=True):
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = False
    sim.slice_step_diagnostics = False
    sim.init_grids()
    if reference_particle:
        sim.beam.ref.set_species("electron").set_kin_energy_MeV(2.0e3)
    if npart > 0:
        sim.add_particles(bunch_charge_C, beam(lambda_pt), npart)
    sim.beam.beam_moments()
    sim.finalize()


moments(npart=0, reference_particle=False)  # a beam a Source loads later
moments(npart=0)  # an empty beam
moments(npart=1000, bunch_charge_C=0.0)  # test particles, carrying no weight
moments(npart=1)  # no emittance in any plane
moments(npart=1000, lambda_pt=0.0)  # no longitudinal emittance

print("no floating-point exception raised")
"""


@pytest.mark.manages_amrex
def test_beam_moments_do_not_raise_fpe(tmp_path):
    """
    Computing the beam moments of a degenerate beam does not raise a floating-point
    exception, i.e. it survives amrex.fpe_trap_invalid / amrex.fpe_trap_zero.
    """
    if not impactx.Config.have_mpi:
        pytest.skip("the stand-alone program pre-initializes MPI via mpi4py")

    program = tmp_path / "fpe_trap_beam_moments.py"
    program.write_text(_FPE_TRAP_PROGRAM)

    result = subprocess.run(
        [sys.executable, str(program)],
        cwd=tmp_path,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    assert "no floating-point exception raised" in result.stdout
