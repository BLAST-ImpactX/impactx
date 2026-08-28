#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
The FFT space-charge solver reuses its Green's function instead of rebuilding it on
every slice step. Reuse only ever happens for a grid the Green's function was built
for, so it must not change the answer at all: these tests compare the reusing solver
against one forced to rebuild every time, and require the results to agree exactly.

Each configuration runs in its own single-threaded subprocess. Charge deposition sums
in a thread-dependent order, so with more than one OpenMP thread even two identical
runs differ in the last digits and an exact comparison would be meaningless. Running
apart also keeps ParmParse settings from leaking between configurations, since the
parameter table outlives an individual simulation.
"""

import json
import math
import os
import subprocess
import sys

import pytest

BUNCH_CHARGE_C = 1.0e-8

# the minimum padding used throughout; the band is set relative to it
PROB_RELATIVE = 1.1


def _simulate(
    rebuild_always,
    max_entries,
    prob_relative_max,
    lattice="constf",
    space_charge="3D",
    n_cell=32,
    nslice=8,
):
    """Run one configuration in-process and return (beam moments, cell size)."""
    import amrex.space3d as amr
    from impactx import ImpactX, distribution, elements

    # set every control on every run: ParmParse outlives a simulation, so a value left
    # over from a previous configuration would silently carry into this one
    pp = amr.ParmParse("ablastr")
    pp.add("igf_rebuild_always", 1 if rebuild_always else 0)
    pp.add("igf_cache_max_entries", max_entries)

    sim = ImpactX()
    # the 2.5D solver projects the charge onto one transverse plane
    sim.n_cell = [n_cell, n_cell, n_cell if space_charge == "3D" else 1]
    if space_charge != "3D":
        sim.blocking_factor_z = [1]
    sim.tiny_profiler = False
    sim.particle_shape = 2
    sim.space_charge = space_charge
    sim.poisson_solver = "fft"
    sim.prob_relative = [PROB_RELATIVE]
    sim.prob_relative_max = prob_relative_max
    sim.slice_step_diagnostics = False
    sim.diagnostics = False
    sim.verbose = 0
    if space_charge == "2p5D":
        # the longitudinal kick gathers the potential itself, not its gradient, so this
        # is the mode in which the additive constant of the 2D Green's function matters
        sim.space_charge_num_longitudinal_bins = 100
        sim.space_charge_apply_longitudinal_kick = True
    sim.init_grids()

    sim.beam.ref.set_species("proton").set_kin_energy_MeV(2.0e3)
    distr = distribution.Waterbag(
        lambdaX=1.2154443728379865788e-3,
        lambdaY=1.2154443728379865788e-3,
        lambdaT=4.0956844276541331005e-4,
        lambdaPx=8.2274435782286157175e-4,
        lambdaPy=8.2274435782286157175e-4,
        lambdaPt=2.4415943602685364584e-3,
    )
    sim.add_particles(BUNCH_CHARGE_C, distr, 10000)
    if lattice == "accel":
        # A uniform accelerating section. ez is normalized, q*Ez/(m*c^2) in 1/m, and 0.4
        # over 2 m raises gamma by about 25%. Both the reference energy and, through
        # adiabatic damping, the beam size then move fast enough to land on a different
        # mesh at practically every step, so this exercises invalidation rather than reuse.
        sim.lattice.extend(
            [elements.ChrAcc(name="acc1", ds=2.0, ez=0.4, bz=0.0, nslice=nslice)]
        )
    else:
        sim.lattice.extend(
            [elements.ConstF(name="cf1", ds=2.0, kx=1.0, ky=1.0, kt=1.0, nslice=nslice)]
        )
    sim.track_particles()

    moments = {k: float(v) for k, v in dict(sim.beam.beam_moments()).items()}
    cell_size = [float(c) for c in sim.Geom(lev=0).data().CellSize()]
    sim.finalize()
    return moments, cell_size


def _run(
    rebuild_always=False,
    max_entries=8,
    prob_relative_max=1.21,
    lattice="constf",
    space_charge="3D",
):
    """Run one configuration in a single-threaded subprocess."""
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    args = [
        sys.executable,
        os.path.abspath(__file__),
        "--emit",
        "1" if rebuild_always else "0",
        str(max_entries),
        str(prob_relative_max),
        lattice,
        space_charge,
    ]
    proc = subprocess.run(args, env=env, capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, (
        f"subprocess failed ({proc.returncode})\n"
        f"--- stdout ---\n{proc.stdout[-2000:]}\n--- stderr ---\n{proc.stderr[-2000:]}"
    )
    # the simulation prints as it tracks; the payload is the final line
    payload = proc.stdout.strip().splitlines()[-1]
    result = json.loads(payload)
    return result["moments"], result["cell_size"]


def _mismatches(a, b):
    return {k: (a[k], b[k]) for k in a if a[k] != b[k]}


def test_the_reference_is_reproducible():
    """The comparisons below are only meaningful if a run reproduces itself."""
    first, _ = _run(rebuild_always=True)
    second, _ = _run(rebuild_always=True)
    assert first == second, (
        "two identical single-threaded runs disagree, so the exact comparisons in this "
        f"file cannot mean anything: {_mismatches(first, second)}"
    )


@pytest.mark.parametrize("max_entries", [0, 8])
def test_reuse_does_not_change_the_result(max_entries):
    """Reusing a Green's function must reproduce rebuilding it, bit for bit.

    ``max_entries=0`` keeps only the Green's function currently in use, exercising the
    reuse of an unchanged grid alone. ``max_entries=8`` additionally exercises putting
    Green's functions aside and taking them back out as the mesh moves between allowed lengths.
    """
    reference, _ = _run(rebuild_always=True)
    reused, _ = _run(rebuild_always=False, max_entries=max_entries)

    assert reference == reused, (
        f"reusing the Green's function changed the result (max_entries={max_entries}): "
        f"{_mismatches(reference, reused)}"
    )


def test_exact_fit_matches_rebuild():
    """With no band, the mesh fits the beam and every step gets its own Green's function.

    This is the configuration in which nothing is ever reused, and it must also leave the
    result untouched.
    """
    reference, _ = _run(rebuild_always=True, prob_relative_max=PROB_RELATIVE)
    reused, _ = _run(rebuild_always=False, prob_relative_max=PROB_RELATIVE)

    assert reference == reused, (
        f"reuse changed the result with the mesh fitted exactly: "
        f"{_mismatches(reference, reused)}"
    )


def test_reuse_survives_a_changing_mesh():
    """A mesh that moves every step must invalidate the Green's function every step.

    Under acceleration both the reference energy and the beam size change quickly, so
    nearly every step lands on a different mesh. Reuse that failed to notice would
    quietly apply a stale Green's function, which is a wrong answer rather than a slow
    one, so this is the case where the guard has to be right.
    """
    reference, _ = _run(rebuild_always=True, lattice="accel")
    reused, _ = _run(rebuild_always=False, lattice="accel")

    assert reference == reused, (
        f"reuse changed the result while the reference energy was changing: "
        f"{_mismatches(reference, reused)}"
    )


@pytest.mark.parametrize("space_charge", ["3D", "2D", "2p5D"])
def test_reuse_holds_for_every_space_charge_model(space_charge):
    """Each solver mode must be unaffected by reuse, 2.5D above all.

    The 2.5D longitudinal kick gathers the transverse potential itself rather than its
    gradient. The 2D integrated Green's function is homogeneous only up to an additive
    constant, and a constant in the potential, invisible to a gradient, goes straight
    into the longitudinal momentum. Reuse must therefore not shift it.
    """
    reference, _ = _run(rebuild_always=True, space_charge=space_charge)
    reused, _ = _run(rebuild_always=False, space_charge=space_charge)

    assert reference == reused, (
        f"reuse changed the {space_charge} result: {_mismatches(reference, reused)}"
    )


def test_two_dimensional_modes_never_reuse_across_scales():
    """2D and 2.5D must reuse a Green's function only at the scale it was built for.

    Rescaling a 3D Green's function is exact. The 2D one picks up an additive constant
    under a change of scale, which the 2.5D longitudinal kick would feel, so the solver
    must decline that reuse rather than correct for it. A mesh exactly twice as large is
    the case a 3D solver would happily reuse, so it is the one to check.
    """
    small, cell_small = _run(space_charge="2p5D")
    # the same beam on a mesh a factor of two larger: a 3D solve would reuse across this
    big, cell_big = _run(space_charge="2p5D", prob_relative_max=2.0 * PROB_RELATIVE)

    ratio = cell_big[0] / cell_small[0]
    assert ratio > 1.0, (
        f"the two runs did not end up on different meshes (ratio {ratio!r}), so this "
        f"test cannot see whether reuse crossed scales"
    )
    reference, _ = _run(
        space_charge="2p5D", rebuild_always=True, prob_relative_max=2.0 * PROB_RELATIVE
    )
    assert big == reference, (
        f"reuse changed the 2.5D result on the wider mesh: {_mismatches(reference, big)}"
    )


@pytest.mark.parametrize("prob_relative_max", [1.32, 1.21, 1.1495])
def test_box_lands_on_a_fit_length(prob_relative_max):
    """The transverse mesh must be one of the allowed lengths, 2**(k/m) m.

    This is what makes a mesh recur: two steps whose padded beam sizes fall between the
    same pair of allowed lengths get the identical mesh, and so the identical Green's
    function. The three bands here work out to 4, 8 and 16 lengths per doubling.
    """
    n_cell = 32
    band = prob_relative_max / PROB_RELATIVE
    fit_lengths_per_doubling = math.ceil(1.0 / math.log2(band))

    _, cell_size = _run(prob_relative_max=prob_relative_max)

    # the transverse mesh is rounded up directly; the longitudinal one is rounded in the
    # frame the solver works in and so carries a gamma, which is not checked here
    for d in (0, 1):
        box_width = cell_size[d] * n_cell
        k = math.log2(box_width) * fit_lengths_per_doubling
        assert math.isclose(k, round(k), abs_tol=1e-9), (
            f"mesh width {box_width!r} along direction {d} is not an allowed length of "
            f"2**(k/{fit_lengths_per_doubling}) m: k = {k!r}"
        )


def test_exactly_fitted_box_is_not_a_fit_length():
    """Guard against the check above passing for the wrong reason."""
    n_cell = 32
    _, cell_size = _run(prob_relative_max=PROB_RELATIVE)

    allowed = []
    for d in (0, 1):
        k = math.log2(cell_size[d] * n_cell) * 8
        allowed.append(math.isclose(k, round(k), abs_tol=1e-9))
    assert not all(allowed), (
        "an exactly fitted mesh happened to be an allowed length, so the check above "
        "cannot tell a restricted mesh from a fitted one"
    )


def test_band_below_prob_relative_is_rejected():
    """The mesh cannot be required to be smaller than the padding asks for."""
    from impactx import ImpactX

    sim = ImpactX()
    sim.prob_relative = [PROB_RELATIVE]
    sim.prob_relative_max = 0.5 * PROB_RELATIVE
    sim.space_charge = "3D"
    sim.poisson_solver = "fft"
    with pytest.raises(Exception):
        sim.init_grids()
    sim.finalize()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--emit":
        rebuild_always = sys.argv[2] == "1"
        max_entries = int(sys.argv[3])
        prob_relative_max = float(sys.argv[4])
        lattice = sys.argv[5]
        space_charge = sys.argv[6]
        moments, cell_size = _simulate(
            rebuild_always,
            max_entries,
            prob_relative_max,
            lattice=lattice,
            space_charge=space_charge,
        )
        print(json.dumps({"moments": moments, "cell_size": cell_size}))
    else:
        test_the_reference_is_reproducible()
        test_reuse_does_not_change_the_result(0)
        test_reuse_does_not_change_the_result(8)
        test_exact_fit_matches_rebuild()
