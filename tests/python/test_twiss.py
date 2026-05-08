#!/usr/bin/env python3
#
# Copyright 2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-
"""
Unit tests for the Wolski-6D lattice Twiss diagnostic.

These tests exercise the algorithmic layers in ``impactx.twiss_lattice`` and
the corresponding C++ bindings:

  * numerical properties of the Wolski parameterization matrix (symplecticity,
    block-diagonalization of the one-turn map),
  * a round-trip between uncoupled Courant-Snyder inputs and the Twiss
    read-out,
  * matched Twiss of a periodic FODO ring,
  * graceful handling of a singular dispersion system,
  * consistency between ``sim.lattice.map_trace(ref)`` and the existing
    ``sim.lattice.transfer_map(ref)`` on a multi-slice lattice,
  * analytic transport of Twiss through a drift in line (transfer-line) mode,
  * a sanity check on the fractional tunes of a stable FODO cell.
"""

import warnings

import numpy as np
import pytest

from impactx import Config, ImpactX, RefPart, elements
from impactx.twiss_lattice import (
    a_from_courant_snyder,
    dispersion_from_m,
    tunes_from_m,
    twiss_from_a,
    wolski_a_from_m,
)


def _symplectic_form():
    """Return the 6x6 canonical symplectic form for ordering
    (x, p_x, y, p_y, t, p_t): block-diagonal with blocks [[0, 1], [-1, 0]]."""
    J = np.zeros((6, 6))
    for i in range(3):
        J[2 * i, 2 * i + 1] = 1.0
        J[2 * i + 1, 2 * i] = -1.0
    return J


def _tol():
    """Relative / absolute tolerance tuned to the build precision.

    ImpactX can be built in either double or single precision; tests must
    accept both. The tuples returned here are used as ``(rtol, atol)``.
    """
    if Config.precision == "SINGLE":
        return 1.0e-4, 1.0e-5
    return 1.0e-9, 1.0e-10


def _build_fodo_lattice(k, l_quad, l_drift):
    """Build a symmetric FODO cell of the form QF/2 -- D -- QD -- D -- QF/2.

    Splitting the focusing quadrupole into two halves at the ends makes the
    cell symmetric about its midpoint. Alpha at the ends is then zero, and
    the matched beta functions agree with the classical thin-lens FODO
    formulas in the short-quad limit.

    Parameters
    ----------
    k : float
        Quadrupole strength in 1/m^2 (focusing in x for ``k > 0``).
    l_quad : float
        Length of the full (non-split) quadrupole in meters.
    l_drift : float
        Length of each drift section in meters.
    """
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Quad(name="qf_half", ds=l_quad / 2.0, k=k, nslice=1),
            elements.Drift(name="d1", ds=l_drift, nslice=1),
            elements.Quad(name="qd", ds=l_quad, k=-k, nslice=1),
            elements.Drift(name="d2", ds=l_drift, nslice=1),
            elements.Quad(name="qf_half2", ds=l_quad / 2.0, k=k, nslice=1),
        ]
    )
    return lattice


def test_wolski_a_symplecticity_and_block_diagonalization():
    """The Wolski parameterization matrix A must be symplectic and must
    block-diagonalize the one-turn map M into rotations by 2*pi*tune."""

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    lattice = _build_fodo_lattice(k=0.5, l_quad=0.2, l_drift=1.0)
    M = lattice.transfer_map(ref).to_numpy()

    A, tunes = wolski_a_from_m(M)
    J = _symplectic_form()

    # (1) A is symplectic: A^T J A = J.
    assert np.allclose(A.T @ J @ A, J, rtol=1e-7, atol=1e-9)

    # (2) A^{-1} M A is block-diagonal.
    A_inv = np.linalg.inv(A)
    M_block = A_inv @ M @ A
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            off_block = M_block[2 * i : 2 * i + 2, 2 * j : 2 * j + 2]
            assert np.allclose(off_block, 0.0, atol=1e-7), (
                f"Off-diagonal 2x2 block ({i},{j}) is not zero after "
                f"diagonalization:\n{off_block}"
            )

    # (3) Each rotational mode's diagonal 2x2 block is a rotation by
    # 2*pi*tune. Parabolic modes (for which wolski_a_from_m reports tune 0
    # and substitutes an identity A-block) are not pure rotations and
    # therefore skipped -- M_block retains the original 2x2 shear for them,
    # e.g. a (t, p_t) drift without an RF cavity.
    for k_mode in range(3):
        if tunes[k_mode] == 0.0:
            continue
        block = M_block[2 * k_mode : 2 * k_mode + 2, 2 * k_mode : 2 * k_mode + 2]
        mu = 2.0 * np.pi * tunes[k_mode]
        expected_rotation = np.array(
            [[np.cos(mu), np.sin(mu)], [-np.sin(mu), np.cos(mu)]]
        )
        assert np.allclose(block, expected_rotation, atol=1e-7), (
            f"Mode {k_mode} diagonal 2x2 block is not a rotation by "
            f"2*pi*{tunes[k_mode]:.6f}:\n{block}"
        )


def test_twiss_uncoupled_courant_snyder_roundtrip():
    """Round-trip: build A from uncoupled Courant-Snyder inputs, read
    beta and alpha back via ``twiss_from_a``, and verify that the diagonal
    Mais-Ripken entries reproduce the inputs and the off-diagonal entries
    are zero."""
    beta_x, alpha_x = 3.0, 0.25
    beta_y, alpha_y = 7.0, -0.5
    beta_t, alpha_t = 1.5, 0.0

    A = a_from_courant_snyder(
        beta_x=beta_x,
        alpha_x=alpha_x,
        beta_y=beta_y,
        alpha_y=alpha_y,
        beta_t=beta_t,
        alpha_t=alpha_t,
    )

    tw = twiss_from_a(A)

    # Diagonal (uncoupled) entries must match the inputs.
    assert np.isclose(tw["beta_1x"], beta_x, atol=1e-12)
    assert np.isclose(tw["alpha_1x"], alpha_x, atol=1e-12)
    assert np.isclose(tw["beta_2y"], beta_y, atol=1e-12)
    assert np.isclose(tw["alpha_2y"], alpha_y, atol=1e-12)
    assert np.isclose(tw["beta_3t"], beta_t, atol=1e-12)
    assert np.isclose(tw["alpha_3t"], alpha_t, atol=1e-12)

    # Off-diagonal (coupling) entries are exactly zero by construction.
    native_plane_of_mode = {1: "x", 2: "y", 3: "t"}
    for mode in (1, 2, 3):
        for plane in ("x", "y", "t"):
            if plane != native_plane_of_mode[mode]:
                assert np.isclose(tw[f"beta_{mode}{plane}"], 0.0, atol=1e-12)


def test_matched_ring_twiss_is_periodic():
    """For a periodic FODO ring, the matched Twiss returned by
    ``sim.twiss()`` must be identical at the cell boundary (s=0 and
    s=L), and the total phase advance along the cell must equal
    2*pi*tune as derived from the one-turn map eigenvalues."""
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = "false"
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    sim.lattice.extend(list(_build_fodo_lattice(k=0.5, l_quad=0.2, l_drift=1.0)))

    tw = sim.twiss()

    # Periodic: first and last sampled values of the transverse Twiss must agree.
    for key in ("beta_1x", "beta_2y", "alpha_1x", "alpha_2y"):
        assert np.isclose(tw[key][0], tw[key][-1], rtol=1e-6, atol=1e-9), (
            f"{key} is not periodic across the cell: {tw[key][0]} vs {tw[key][-1]}"
        )

    # Tunes from the one-turn eigendecomposition must match the net phase
    # advance accumulated along the cell (mu_k / (2 pi) == nu_k).
    assert "tunes" in tw
    assert np.isclose(tw["mu_1"][-1] / (2.0 * np.pi), tw["tunes"][0], atol=1e-7)
    assert np.isclose(tw["mu_2"][-1] / (2.0 * np.pi), tw["tunes"][1], atol=1e-7)

    sim.finalize()


def test_parabolic_longitudinal_mode_is_masked_to_nan():
    """A ring without an RF cavity has a longitudinal one-turn 2x2 block
    of the form [[1, c], [0, 1]] (a drift in arrival time). This block is
    defective: both eigenvalues equal one, the Jordan form is non-trivial,
    and no Wolski-type A-block with a finite beta exists. ``compute_twiss``
    must therefore report NaN for ``beta_3t``, ``alpha_3t`` and ``mu_3``
    (and, through the alias, ``beta_t``, ``alpha_t``) along the whole
    lattice, rather than propagating the defective 2x2 shear and returning
    shear-driven numbers that look like valid Twiss values.

    Transverse outputs and dispersion must remain finite.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = "false"
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    sim.lattice.extend(list(_build_fodo_lattice(k=0.5, l_quad=0.2, l_drift=1.0)))

    with pytest.warns(RuntimeWarning, match="not rotational"):
        tw = sim.twiss()

    # Longitudinal Twiss is masked.
    for key in ("beta_3t", "alpha_3t", "mu_3", "beta_t", "alpha_t"):
        assert np.all(np.isnan(tw[key])), (
            f"{key} should be NaN throughout for a parabolic longitudinal "
            f"mode, got {tw[key]}"
        )
    # The longitudinal tune is reported as exactly zero.
    assert tw["tunes"][2] == 0.0

    # Transverse Twiss and transverse phase advance must remain finite.
    for key in (
        "beta_1x",
        "alpha_1x",
        "mu_1",
        "beta_2y",
        "alpha_2y",
        "mu_2",
        "beta_x",
        "beta_y",
        "alpha_x",
        "alpha_y",
    ):
        assert np.all(np.isfinite(tw[key])), f"{key} must remain finite, got {tw[key]}"

    # Transverse-in-transverse coupling entries are finite too; the only
    # masked entries are the three columns tied to the parabolic mode.
    for key in ("beta_1y", "beta_2x", "beta_1t", "beta_2t"):
        assert np.all(np.isfinite(tw[key])), f"{key} must remain finite, got {tw[key]}"

    # Dispersion is unaffected by the longitudinal mask.
    for key in ("disp_x", "disp_px", "disp_y", "disp_py"):
        assert np.all(np.isfinite(tw[key])), f"{key} must remain finite, got {tw[key]}"

    sim.finalize()


def test_dispersion_is_zero_for_pure_quad_drift_fodo():
    """A FODO of only quadrupoles and drifts has no dispersion at the
    closed orbit: the coupling column M[0:4, 5] is zero, so the periodic
    equation ``(I - M_4x4) D = 0`` has D = 0 as the unique solution (since
    the transverse tune is strictly between integer values, (I - M_4x4)
    is non-singular for a stable FODO)."""
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    lattice = _build_fodo_lattice(k=0.5, l_quad=0.2, l_drift=1.0)
    M = lattice.transfer_map(ref).to_numpy()

    D = dispersion_from_m(M)
    assert np.allclose(D, 0.0, atol=1e-12)


def test_map_trace_and_transfer_map_agree():
    """Regression test for the slice-accumulation convention.

    For a mixed multi-slice lattice (Sbend, Drift, Quad, Drift, each with
    ``nslice=2``), both public APIs

        sim.lattice.map_trace(ref)[-1]["M"]  -- last entry of the per-element trace
        sim.lattice.transfer_map(ref)        -- end-to-end linear transfer map

    must return the same 6x6 matrix, and that matrix must agree with the
    reference values locked in by ``test_lattice_optics.py``. Since
    ``transfer_map`` is implemented as a thin wrapper over
    ``diagnostics::linear_map`` (shared with ``map_trace``), any future
    change to the per-slice ordering will show up here.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = "false"
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    # Same lattice as in tests/python/test_lattice_optics.py::test_lattice_linear_map
    # so the hard-coded expected R below remains a single point of truth.
    sim.lattice.extend(
        [
            elements.Sbend(name="bend1", ds=1.0, rc=10.0, nslice=2),
            elements.Drift(name="drift1", ds=2.0, nslice=2),
            elements.Quad(name="quad1", ds=0.5, k=1.0, nslice=2),
            elements.Drift(name="drift2", ds=1.0, nslice=2),
        ]
    )

    R_expected = np.array(
        [
            [3.74670546e-01, 2.54005827e00, 0, 0, 0, -2.34864805e-01],
            [-4.76219076e-01, -5.59489407e-01, 0, 0, 0, 3.20646253e-02],
            [0, 0, 1.64872127e00, 6.59488508e00, 0, 0],
            [0, 0, 5.21095305e-01, 2.69091188e00, 0, 0],
            [9.98334297e-02, 4.99583537e-02, 0, 0, 1.00000000e00, -1.66466013e-03],
            [0, 0, 0, 0, 0, 1.00000000e00],
        ]
    )

    if Config.precision == "SINGLE":
        atol, rtol = 1.0e-7, 5.0e-5
    else:
        atol, rtol = 0.0, 1.0e-8

    trace = sim.lattice.map_trace(sim.beam.ref)
    # One leading identity entry, one entry per lattice element.
    assert len(trace) == len(sim.lattice) + 1

    M_trace_final = trace[-1]["M"].to_numpy()
    M_transfer_map = sim.lattice.transfer_map(sim.beam.ref).to_numpy()

    # The two APIs must agree on the end-to-end map.
    assert np.allclose(M_trace_final, M_transfer_map, rtol=rtol, atol=max(atol, 1e-8))

    # And the agreed-upon value must match the reference from test_lattice_optics.py.
    assert np.allclose(M_transfer_map, R_expected, rtol=rtol, atol=atol)

    # The leading entry of the trace must be the identity map at the start.
    assert np.allclose(trace[0]["M"].to_numpy(), np.eye(6), atol=1e-15)

    sim.finalize()


def test_lattice_diagnostics_do_not_mutate_caller_reference():
    """``sim.lattice.transfer_map(ref)``, ``sim.lattice.map_trace(ref)``,
    and the high-level ``sim.twiss()`` wrapper must not advance the
    caller's reference particle in place.

    Users commonly inspect linear optics before deciding to run a tracking
    simulation afterward; an accidental in-place advance would silently
    shift the starting ``s`` and drift coordinates of the reference and
    change subsequent tracking results.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = "false"
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    sim.lattice.extend(
        [
            elements.Drift(name="d", ds=1.0, nslice=2),
            elements.Quad(name="q", ds=0.5, k=0.4, nslice=2),
        ]
    )

    # Snapshot the reference-particle fields whose mutation would be
    # physically observable (``s``, ``x``, ``px``, ``z``, ``t``, ``pt``).
    before = (ref.s, ref.x, ref.px, ref.z, ref.t, ref.pt)

    sim.lattice.transfer_map(sim.beam.ref)
    assert (ref.s, ref.x, ref.px, ref.z, ref.t, ref.pt) == before

    sim.lattice.map_trace(sim.beam.ref)
    assert (ref.s, ref.x, ref.px, ref.z, ref.t, ref.pt) == before

    sim.twiss()
    assert (ref.s, ref.x, ref.px, ref.z, ref.t, ref.pt) == before

    sim.finalize()


def test_transfer_map_handles_sedge_reset_with_soft_element():
    """Regression test for the ``ref.sedge`` reset at each element entry.

    The ``SoftSolenoid``, ``SoftQuadrupole``, and ``RFCavity`` elements
    compute their local integration coordinate as
    ``zin = ref.s - ref.sedge``. For this to be correct at the second
    (and any subsequent) edge-sensitive element in a lattice,
    ``transfer_map`` must reset ``ref.sedge = ref.s`` at each element
    entry, matching the tracking loop in ``src/tracking/envelope.cpp``.

    Before this fix the default starting ``ref.sedge = 0`` was never
    refreshed, so the second SoftSolenoid would see
    ``zin = L_sol_1 + L_sol_2`` (way outside its own field profile)
    instead of ``zin = 0`` at its entry.

    For a static-field axially-varying element such as ``SoftSolenoid``
    the per-element R-matrix only depends on the reference momentum
    (invariant through the element) and the local ``zin``. Therefore
    concatenating two copies of the same element must give the square
    of the single-element R-matrix:

        transfer_map([A, A]) == transfer_map([A]) @ transfer_map([A])

    If ``ref.sedge`` were *not* reset at the second A's entry, the R of
    the second A would be wrong and the equality would fail.
    """
    # A small one-lobe-cosine field profile. The numerics of the
    # profile are unimportant here; we only need the element to actually
    # read ref.sedge in its per-slice push.
    cos_coefficients = [1.0, 0.2, 0.05]
    sin_coefficients = [0.0, 0.0, 0.0]

    def make_lattice(n_copies):
        lattice = elements.KnownElementsList()
        lattice.extend(
            [
                elements.SoftSolenoid(
                    name=f"sol{i}",
                    ds=0.3,
                    bscale=1.0,
                    cos_coefficients=cos_coefficients,
                    sin_coefficients=sin_coefficients,
                    unit=0,
                    mapsteps=20,
                    nslice=2,
                )
                for i in range(n_copies)
            ]
        )
        return lattice

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    M_one = make_lattice(1).transfer_map(ref).to_numpy()
    M_two = make_lattice(2).transfer_map(ref).to_numpy()

    # Concatenation identity: (static-field, momentum-preserving element) -> R^2.
    assert np.allclose(M_two, M_one @ M_one, atol=1e-10), (
        "transfer_map([A, A]) != transfer_map([A])^2 for a static "
        "SoftSolenoid. The second SoftSolenoid likely saw a stale "
        "ref.sedge from before the first element, producing a wrong "
        "per-slice R."
    )


def test_transfer_map_fallback_kwarg_does_not_warn():
    """When the user opts in to ``fallback_identity_map=True``, the
    identity substitution must be silent (no ``RuntimeWarning`` nor
    ABLASTR warning). The opt-in is an explicit acknowledgment that
    non-linear elements are present, so warning on every call would be
    noise. Without the opt-in (default), the call must raise."""
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(name="d", ds=0.5, nslice=1),
            # NonlinearLens has no linear transport map by design.
            elements.NonlinearLens(knll=1.0e-3, cnll=0.1, name="nll"),
            elements.Drift(name="d2", ds=0.5, nslice=1),
        ]
    )

    # Default behavior: must raise because a NonlinearLens has no linear map.
    with pytest.raises(RuntimeError, match="Undefined linear transport map"):
        lattice.transfer_map(ref)

    # Explicit opt-in: must succeed silently (no warnings emitted).
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning becomes an exception
        M = lattice.transfer_map(ref, fallback_identity_map=True).to_numpy()

    # With NonlinearLens treated as identity, the lattice reduces to two
    # drifts of total length 1.0 m.
    expected = np.eye(6)
    expected[0, 1] = 1.0
    expected[2, 3] = 1.0
    expected[4, 5] = 1.0 / (ref.pt**2 - 1.0)  # drift longitudinal term
    assert np.allclose(M[0:4, 0:4], expected[0:4, 0:4], atol=1e-10)


def test_dispersion_raises_on_singular_lattice():
    """When ``(I - M_4x4)`` is singular (integer transverse tune or an
    unstable lattice), ``dispersion_from_m`` must raise a
    ``LinAlgError`` with a message that identifies the cause, instead of
    silently returning numerically meaningless values.

    A single Marker element is used here: its 4x4 transverse block is
    the identity, so ``I - M_4x4`` is the zero matrix -- the degenerate
    limit of an integer tune.
    """
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    lattice = elements.KnownElementsList()
    lattice.extend([elements.Marker(name="m0")])

    M = lattice.transfer_map(ref).to_numpy()
    with pytest.raises(np.linalg.LinAlgError, match="singular|integer"):
        dispersion_from_m(M)


def test_line_twiss_propagates_drift_analytically():
    """Propagating Twiss through a drift of length L with initial
    ``(beta0, alpha0 = 0)`` must reproduce the textbook formulas

        beta(s)  = beta0 + s^2 / beta0,
        alpha(s) = -s / beta0.

    This validates the line-mode (non-periodic) branch of
    ``compute_twiss`` against a closed-form result.
    """
    sim = ImpactX()
    sim.particle_shape = 2
    sim.space_charge = "false"
    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    L = 2.0
    sim.lattice.extend([elements.Drift(name="long_drift", ds=L, nslice=1)])

    beta0, alpha0 = 4.0, 0.0

    tw = sim.twiss(
        init={
            "beta_x": beta0,
            "alpha_x": alpha0,
            "beta_y": beta0,
            "alpha_y": alpha0,
        }
    )

    expected_beta = beta0 + L * L / beta0
    expected_alpha = -L / beta0

    assert np.isclose(tw["beta_1x"][-1], expected_beta, rtol=1e-7, atol=1e-9)
    assert np.isclose(tw["alpha_1x"][-1], expected_alpha, rtol=1e-7, atol=1e-9)
    assert np.isclose(tw["beta_2y"][-1], expected_beta, rtol=1e-7, atol=1e-9)
    assert np.isclose(tw["alpha_2y"][-1], expected_alpha, rtol=1e-7, atol=1e-9)

    sim.finalize()


def test_fodo_tunes_in_stable_range():
    """For a symmetric FODO cell chosen in the stable region of the
    Courant-Snyder stability diagram, the two transverse fractional
    tunes must lie strictly in the open interval (0, 1/2)."""
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    lattice = _build_fodo_lattice(k=0.8, l_quad=0.2, l_drift=0.8)
    M = lattice.transfer_map(ref).to_numpy()

    tunes = tunes_from_m(M)
    for nu in tunes[:2]:
        assert 0.0 < abs(nu) < 0.5, (
            f"transverse tune {nu} is outside the stable range (0, 1/2)"
        )
