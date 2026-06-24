#!/usr/bin/env python3
#
# Copyright 2022-2024 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, RefPart, elements


def test_lattice_linear_map():
    """Calculate the linear transfer map of a lattice."""

    # Create reference particle
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    # Create a valid test lattice (all elements define a linear transfer map)
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Sbend(name="bend1", ds=1.0, rc=10.0, nslice=2),
            elements.Drift(name="drift1", ds=2.0, nslice=2),
            elements.Quad(name="quad1", ds=0.5, k=1.0, nslice=2),
            elements.Drift(name="drift2", ds=1.0, nslice=2),
        ]
    )

    # Expected result (matrix multiplication)
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
        atol = 1.0e-7
        rtol = 5.0e-5
    else:
        atol = 0.0
        rtol = 1.0e-8

    # Calculate Linear Transfer Map
    R = lattice.transfer_map(ref)
    assert np.allclose(R.to_numpy(), R_expected, rtol=rtol, atol=atol)

    # Check unexpected/unsupported options
    with pytest.raises(RuntimeError):
        lattice.transfer_map(ref, order="invalid")

    # Create a lattice with an element that does not define a linear transfer
    # map. Programmable has no built-in linear map; users may opt in via
    # its transport_map_hook attribute (see test_programmable_linear_hook).
    lattice.append(elements.Programmable())

    # Ensure that the calculation asserts
    with pytest.raises(RuntimeError):
        lattice.transfer_map(ref)

    # Now the user explicitly assumes that undefined maps are identity maps
    R = lattice.transfer_map(ref, fallback_identity_map=True)
    assert np.allclose(R.to_numpy(), R_expected, rtol=rtol, atol=atol)


def _symplectic_form_6x6():
    """Standard 6x6 symplectic form for the (x, px, y, py, t, pt) ordering."""
    J = np.zeros((6, 6))
    for i in range(3):
        J[2 * i, 2 * i + 1] = 1.0
        J[2 * i + 1, 2 * i] = -1.0
    return J


def test_mixed_lattice_is_symplectic():
    """Compose several newly-implemented element maps and check that the
    end-to-end 6x6 linear map is symplectic.

    This is a single end-to-end check that every new element's map is
    individually symplectic and that they compose correctly. It avoids
    re-encoding any of the per-element C++ formulas.
    """
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(250.0)

    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(ds=0.4),
            elements.Multipole(multipole=2, K_normal=1.5, K_skew=0.3),
            elements.Drift(ds=0.25),
            # Identity contribution; also exercises the aperture element
            # in a linear-optics call.
            elements.PolygonAperture(
                vertices_x=[-0.02, 0.02, 0.02, -0.02, -0.02],
                vertices_y=[-0.02, -0.02, 0.02, 0.02, -0.02],
            ),
            elements.Drift(ds=0.3),
            elements.NonlinearLens(knll=0.05, cnll=0.25),
            elements.Drift(ds=0.3),
            elements.TaperedPL(k=1.2, taper=0.0, unit=0),
            elements.Drift(ds=0.3),
            elements.ThinDipole(theta=1.5, rc=8.0),
            elements.Drift(ds=0.4),
            elements.PRot(phi_in=2.0, phi_out=-1.0),
            elements.Drift(ds=0.4),
        ]
    )

    R = lattice.transfer_map(ref).to_numpy()
    J = _symplectic_form_6x6()

    # Symplectic condition R^T J R = J; ChrAcc is not in this lattice so
    # the overall composition is expected to be symplectic.
    violation = R.T @ J @ R - J
    tol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-10
    assert np.max(np.abs(violation)) < tol, (
        f"R^T J R - J has max abs entry {np.max(np.abs(violation))}, exceeds tol {tol}."
    )


def test_nonlinear_lens_vanishes_in_linear_limit():
    """At knll -> 0, the NonlinearLens linear map must reduce to the
    identity. This guards against a stray O(knll^2) term leaking into
    the linearization."""
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(100.0)

    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.Drift(ds=0.5),
            elements.NonlinearLens(knll=0.0, cnll=0.4),
            elements.Drift(ds=0.5),
        ]
    )
    R_nl = lattice.transfer_map(ref).to_numpy()

    lattice_ref = elements.KnownElementsList()
    lattice_ref.extend([elements.Drift(ds=1.0)])
    R_drift = lattice_ref.transfer_map(ref).to_numpy()

    tol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-12
    assert np.max(np.abs(R_nl - R_drift)) < tol


def test_chracc_adiabatic_damping_invariant():
    """Uniform-Ez acceleration without a solenoid (Bz = 0) must preserve
    the adiabatic-damping invariant <x^2> * (bg)^2 in the diagonal 1,1
    entry of the transfer map, up to the unit-conversion factor.

    In ImpactX's static phase-space convention (px normalized by the
    local reference momentum), a pure acceleration multiplies every
    px-column entry by bgi/bgf. The (x, x) entry is unaffected and the
    (x, px) entry is bgi/bgf times the drift-like length integral. The
    physical invariant is then (x, x_physical) = x * bgf^0 = x; and
    (sigma_x * sigma_px) scales as (bgi/bgf) as it should for adiabatic
    damping of the geometric emittance.

    Concretely: R(1, 1) = 1, R(2, 2) = bgi/bgf (no solenoid, no
    focusing), R(2, 1) = 0 and R(6, 6) = bgi/bgf. This isolates the
    static<->dynamic unit change without the solenoid rotation and is
    the strongest analytic cross-check available for ChrAcc's linear map.
    """
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(50.0)

    # Small Ez so that theta = alpha/ez * log(...) = 0 when Bz = 0 anyway
    # (alpha = Bz/2 = 0). We still exercise the static<->dynamic rescaling.
    ez = 0.5  # normalized, 1/m
    ds = 0.75  # m
    lattice = elements.KnownElementsList()
    lattice.extend(
        [
            elements.ChrAcc(ds=ds, ez=ez, bz=0.0, nslice=4),
        ]
    )
    R = lattice.transfer_map(ref).to_numpy()

    # Recompute bgi (initial, before the lattice) and bgf (final, after
    # accelerating through ds) from the reference-particle convention.
    pti = ref.pt  # initial pt (negative: pt = -gamma)
    ptf = pti - ez * ds
    bgi = np.sqrt(pti**2 - 1.0)
    bgf = np.sqrt(ptf**2 - 1.0)

    tol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-10

    # (x, x) must stay 1 (transverse position is unchanged by pure Ez).
    assert abs(R[0, 0] - 1.0) < tol
    # No transverse focusing from Bz = 0 means R(2, 1) = R(4, 3) = 0.
    assert abs(R[1, 0]) < tol
    assert abs(R[3, 2]) < tol
    # Adiabatic damping: px_out_static = (bgi/bgf) * px_in_static.
    assert abs(R[1, 1] - bgi / bgf) < tol * (bgi / bgf)
    assert abs(R[3, 3] - bgi / bgf) < tol * (bgi / bgf)
    # pt has the same bgi/bgf scaling in the static unit convention.
    assert abs(R[5, 5] - bgi / bgf) < tol * (bgi / bgf)


def test_programmable_linear_hook():
    """When Programmable.transport_map_hook is set, it must be used by
    transfer_map. When it is unset, transfer_map must raise."""
    from impactx import Map6x6

    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(100.0)

    # A user-supplied linear map: a thin quad kick with strength k_user.
    k_user = 0.7
    R_user = np.eye(6)
    R_user[1, 0] = -k_user
    R_user[3, 2] = -k_user

    def hook(_refpart):
        out = Map6x6()
        for i in range(6):
            for j in range(6):
                out[i + 1, j + 1] = float(R_user[i, j])
        return out

    # Without hook: must raise. (extend() copies the element into the
    # lattice's variant, so any subsequent modification of the local
    # Python handle does not propagate into the lattice.)
    lattice_no_hook = elements.KnownElementsList()
    lattice_no_hook.extend(
        [elements.Drift(ds=0.1), elements.Programmable(), elements.Drift(ds=0.1)]
    )
    with pytest.raises(RuntimeError):
        lattice_no_hook.transfer_map(ref)

    # With hook set before insertion: drift(0.1) * R_user * drift(0.1).
    prog = elements.Programmable()
    prog.transport_map_hook = hook
    lattice = elements.KnownElementsList()
    lattice.extend([elements.Drift(ds=0.1), prog, elements.Drift(ds=0.1)])
    R = lattice.transfer_map(ref).to_numpy()

    pt = ref.pt
    bg2 = pt * pt - 1.0
    D = np.eye(6)
    D[0, 1] = 0.1
    D[2, 3] = 0.1
    D[4, 5] = 0.1 / bg2
    R_expected = D @ R_user @ D

    tol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-10
    assert np.max(np.abs(R - R_expected)) < tol


def test_analytic_element_entries():
    """Analytic spot checks for the signs and special cases of each new
    element map. This is deliberately compact: one lattice per element,
    each with a closed-form expected value that's easy to eyeball.
    """
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(500.0)
    pt = ref.pt
    beta = np.sqrt(1.0 - 1.0 / pt**2)

    tol = 1.0e-4 if Config.precision == "SINGLE" else 1.0e-10

    def R_of(*els):
        lat = elements.KnownElementsList()
        lat.extend(list(els))
        return lat.transfer_map(ref).to_numpy()

    # --- NonlinearLens with a nonzero knll: thin quad with opposite-sign
    # focusing in the two transverse planes ---
    knll, cnll = 0.08, 0.3
    R_nl = R_of(elements.NonlinearLens(knll=knll, cnll=cnll))
    k_lin = 2.0 * knll / cnll**2
    assert abs(R_nl[1, 0] - (-k_lin)) < tol, "NonlinearLens R(2,1)"
    assert abs(R_nl[3, 2] - (+k_lin)) < tol, "NonlinearLens R(4,3)"

    # --- ThinDipole signs: R(2,1)=-theta/rc, R(2,6)=-theta/beta,
    # R(5,1)=+theta/beta ---
    theta_deg, rc = 1.5, 6.0
    theta = np.deg2rad(theta_deg)
    R_td = R_of(elements.ThinDipole(theta=theta_deg, rc=rc))
    assert abs(R_td[1, 0] - (-theta / rc)) < tol, "ThinDipole R(2,1)"
    assert abs(R_td[1, 5] - (-theta / beta)) < tol, "ThinDipole R(2,6)"
    assert abs(R_td[4, 0] - (+theta / beta)) < tol, "ThinDipole R(5,1)"

    # --- PRot inverse: PRot(a, b) then PRot(b, a) must give identity ---
    R_prot = R_of(
        elements.PRot(phi_in=3.0, phi_out=-1.5),
        elements.PRot(phi_in=-1.5, phi_out=3.0),
    )
    assert np.max(np.abs(R_prot - np.eye(6))) < tol, "PRot composes to identity"

    # --- Multipole identity cases (m==1 dipole kick; m>=3) ---
    assert (
        np.max(
            np.abs(
                R_of(elements.Multipole(multipole=1, K_normal=2.7, K_skew=1.3))
                - np.eye(6)
            )
        )
        < tol
    )
    assert (
        np.max(
            np.abs(
                R_of(elements.Multipole(multipole=3, K_normal=7.0, K_skew=0.5))
                - np.eye(6)
            )
        )
        < tol
    )
    assert (
        np.max(
            np.abs(
                R_of(elements.Multipole(multipole=5, K_normal=9.0, K_skew=0.2))
                - np.eye(6)
            )
        )
        < tol
    )
    # and m==2 matches the thin quad
    Kn, Ks = 1.1, -0.4
    R_m2 = R_of(elements.Multipole(multipole=2, K_normal=Kn, K_skew=Ks))
    assert abs(R_m2[1, 0] - (-Kn)) < tol
    assert abs(R_m2[1, 2] - (+Ks)) < tol
    assert abs(R_m2[3, 0] - (+Ks)) < tol
    assert abs(R_m2[3, 2] - (+Kn)) < tol

    # --- TaperedPL: taper is quadratic -> must not enter the linear map ---
    k_val = 1.3
    R_taper0 = R_of(elements.TaperedPL(k=k_val, taper=0.0, unit=0))
    R_taper1 = R_of(elements.TaperedPL(k=k_val, taper=2.5, unit=0))
    assert np.max(np.abs(R_taper0 - R_taper1)) < tol, (
        "taper must not enter the linear map"
    )
    assert abs(R_taper0[1, 0] - (-k_val)) < tol
    assert abs(R_taper0[3, 2] - (-k_val)) < tol
    # unit=1: g = k / rigidity_Tm -> scale by rigidity
    rigidity = ref.rigidity_Tm
    k_T = 0.5  # in T for unit=1
    R_unit1 = R_of(elements.TaperedPL(k=k_T, taper=0.0, unit=1))
    assert abs(R_unit1[1, 0] - (-k_T / rigidity)) < tol
    assert abs(R_unit1[3, 2] - (-k_T / rigidity)) < tol


def test_chracc_solenoid_limit_bz_nonzero():
    """ChrAcc with nonzero Bz reduces to a solenoid in the weak-Ez limit.

    Unit convention difference between the two elements:
      ChrAcc's  bz = q*Bz / (m*c)  (1/m, absolute normalization)
      Sol's     ks = q*Bz / p_ref = bz / (beta*gamma)_ref  (1/m, static norm.)
    so in the ez -> 0 limit (no acceleration), ChrAcc(bz=B) must match
    Sol(ks = B / betagamma_ref), not Sol(ks=B). We check the entire
    transverse 4x4 block as one analytic invariant.
    """
    ref = RefPart()
    ref.set_species("electron").set_kin_energy_MeV(500.0)

    ds = 0.2
    bz = 3.0  # ChrAcc solenoid strength in 1/m (per m*c)
    ez = 1.0e-5  # small Ez, so ChrAcc -> Sol

    lat_chracc = elements.KnownElementsList()
    lat_chracc.extend([elements.ChrAcc(ds=ds, ez=ez, bz=bz, nslice=1)])
    R_c = lat_chracc.transfer_map(ref).to_numpy()

    # Convert ChrAcc's bz (per m*c) into Sol's ks (per rigidity) using the
    # reference betagamma at the entrance of the slice.
    pt = ref.pt
    bg0 = np.sqrt(pt * pt - 1.0)
    ks_equiv = bz / bg0

    lat_sol = elements.KnownElementsList()
    lat_sol.extend([elements.Sol(ds=ds, ks=ks_equiv, nslice=1)])
    R_s = lat_sol.transfer_map(ref).to_numpy()

    # In the ez -> 0 limit, ChrAcc's 4x4 transverse block must match Sol's.
    # Residual O(ez*ds) differences (adiabatic damping, reference-bg drift)
    # are well below the floor set here.
    diff = np.max(np.abs(R_c[:4, :4] - R_s[:4, :4]))
    tol = 5.0e-4 if Config.precision == "SINGLE" else 1.0e-6
    assert diff < tol, (
        f"ChrAcc(ez->0, bz) does not reduce to Sol(ks=bz/betagamma); "
        f"max 4x4 diff = {diff}"
    )
