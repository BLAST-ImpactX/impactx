#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# This set of tests verifies that the particle map of each element is
# symplectic: its Jacobian M satisfies M^T J M = J. Where an element changes
# the reference energy, the momenta are normalized to a different reference
# momentum after the element, and M^T J M = lambda J, with lambda the ratio of
# the reference momenta before and after the element.
#
# -*- coding: utf-8 -*-

from math import factorial

import numpy as np
import pytest

from impactx import Config, ImpactX, Map6x6, elements

# element config
nslice = 1
mapsteps = 4

# no apertures: the particles are far off axis
ALIGNMENT_KWARGS = dict(dx=1.0e-3, dy=-2.0e-3, rotation=0.5)

# canonical pairs (x, px), (y, py), (t, pt)
PHASE_COLS = [
    "position_x",
    "momentum_x",
    "position_y",
    "momentum_y",
    "position_t",
    "momentum_t",
]
J = np.kron(np.eye(3), np.array([[0.0, 1.0], [-1.0, 0.0]]))

# a particle far off axis, where nonlinear fields are strong: x, px, y, py, t, pt
CENTER = np.array([0.05, 0.01, -0.03, 0.005, 0.0, 0.0])

# finite-difference step of the Jacobian, and tolerance of M^T J M - lambda J
if Config.precision == "SINGLE":
    step = 1.0e-3
    atol = 3.0e-4
else:
    step = 1.0e-4
    atol = 2.0e-9


# Elements are tested with and without spin only if their spin push integrates
# the orbit in its own code path; otherwise, it calls the same orbit push.
@pytest.fixture(scope="function")
def sim(request):
    sim = ImpactX()
    sim.space_charge = False
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    sim.spin = getattr(request, "param", False)  # default to False if not parametrized
    sim.init_grids()

    sim.beam.ref.set_species("electron").set_kin_energy_MeV(1.0e3)

    yield sim

    sim.finalize()


def jacobian(el, sim, center=CENTER):
    """Jacobian of the particle map of an element, by fourth-order central differences

    @return the Jacobian and the ratio of the reference momenta before and after
    """
    offsets = [2.0, 1.0, -1.0, -2.0]
    phase_space = np.vstack([center + o * step * np.eye(6) for o in offsets])
    x, px, y, py, t, pt = phase_space.T
    spin = {}
    if sim.spin:
        s = np.full_like(x, 1.0 / np.sqrt(3.0))
        spin = dict(sx=s, sy=s, sz=s)

    beam = sim.beam
    beam.add_n_particles(x, y, t, px, py, pt, beam.ref.qm_ratio_SI, 1.0e-12, **spin)
    beta_gamma_in = beam.ref.beta_gamma
    el.push(beam)
    beta_gamma_out = beam.ref.beta_gamma

    # lost particles stay in the container, flagged invalid
    n_valid = beam.total_number_of_particles(only_valid=True, only_local=True)
    assert n_valid == len(phase_space), "particles were lost"

    df = beam.to_df(local=True).reset_index(drop=True).sort_values("idcpu")
    plus2, plus1, minus1, minus2 = np.split(
        df[PHASE_COLS].to_numpy(dtype=np.float64), len(offsets)
    )
    M = (-plus2 + 8.0 * plus1 - 8.0 * minus1 + minus2).T / (12.0 * step)
    return M, beta_gamma_in / beta_gamma_out


def check_symplectic(el, sim, center=CENTER, atol=atol):
    """Assert that the particle map of an element is symplectic at a point"""
    M, ratio = jacobian(el, sim, center)
    np.testing.assert_allclose(M.T @ J @ M, ratio * J, rtol=0.0, atol=atol)


# =============================================================================
# Thick elements
# =============================================================================


def test_CFbend(sim):
    check_symplectic(
        elements.CFbend(
            ds=0.5, rc=7.613657587094493, k=-7.057403, nslice=nslice, **ALIGNMENT_KWARGS
        ),
        sim,
    )


def test_ChrDrift(sim):
    check_symplectic(elements.ChrDrift(ds=0.25, nslice=nslice, **ALIGNMENT_KWARGS), sim)


@pytest.mark.parametrize(
    ("unit", "k"),
    [(0, 2.98636067687944129), (0, -2.98636067687944129), (1, 5.0)],
    ids=["madx", "madx-defocusing", "si"],
)
def test_ChrPlasmaLens(sim, unit, k):
    check_symplectic(
        elements.ChrPlasmaLens(
            ds=0.331817852986604588, k=k, unit=unit, nslice=nslice, **ALIGNMENT_KWARGS
        ),
        sim,
    )


@pytest.mark.parametrize(("unit", "k"), [(0, 1.0), (1, 3.5)], ids=["madx", "marylie"])
def test_ChrQuad(sim, unit, k):
    check_symplectic(
        elements.ChrQuad(ds=1.0, k=k, unit=unit, nslice=nslice, **ALIGNMENT_KWARGS),
        sim,
    )


def test_ChrAcc(sim):
    check_symplectic(
        elements.ChrAcc(
            ds=1.8,
            ez=10871.950994502130424,
            bz=250.0,
            nslice=nslice,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


def test_ConstF(sim):
    # M[pt, t] = -kt (beta gamma)^2 sin(kt ds) is of order one for this kt, which
    # keeps the energies of the finite-difference stencil physical
    ds = 2.0
    kt = 1.0 / (sim.beam.ref.beta_gamma * np.sqrt(ds))
    check_symplectic(
        elements.ConstF(
            ds=ds, kx=1.0, ky=1.0, kt=kt, nslice=nslice, **ALIGNMENT_KWARGS
        ),
        sim,
    )


def test_Drift(sim):
    check_symplectic(elements.Drift(ds=0.25, nslice=nslice, **ALIGNMENT_KWARGS), sim)


def test_ExactDrift(sim):
    check_symplectic(
        elements.ExactDrift(ds=0.25, nslice=nslice, **ALIGNMENT_KWARGS), sim
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
@pytest.mark.parametrize(
    ("unit", "k_normal", "k_skew"),
    [
        (0, [0.0, 1.0, 0.15], [0.0, 0.35, -0.05]),
        (1, [0.0, 3.5, 0.45], [0.0, 0.8, -0.12]),
    ],
    ids=["madx", "si"],
)
def test_ExactMultipole(sim, unit, k_normal, k_skew):
    check_symplectic(
        elements.ExactMultipole(
            ds=1.0,
            k_normal=k_normal,
            k_skew=k_skew,
            unit=unit,
            mapsteps=mapsteps,
            nslice=nslice,
            int_order=4,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
@pytest.mark.parametrize(
    ("unit", "k_normal", "k_skew"),
    [
        (0, [0.12, 0.015, -0.002], [0.0, 0.01, 0.001]),
        (1, [0.4, 0.05, -0.008], [0.0, 0.03, 0.004]),
    ],
    ids=["madx", "si"],
)
def test_ExactCFbend(sim, unit, k_normal, k_skew):
    tol = atol
    if Config.precision == "SINGLE":
        # FIXME: the bend map cancels terms of order rc in single precision,
        # M^T J M - J is up to 9e-3 with spin and 4e-3 without, which hides
        # smaller defects here; test_ExactCFbend_multipoles checks each multipole
        # order to 1e-3. With #1774 and #1776, it is up to 2e-5 and the default
        # atol holds.
        # https://github.com/BLAST-ImpactX/impactx/issues/1771
        tol = 2.0e-2
    check_symplectic(
        elements.ExactCFbend(
            ds=1.0,
            k_normal=k_normal,
            k_skew=k_skew,
            unit=unit,
            int_order=4,
            mapsteps=mapsteps,
            nslice=nslice,
            **ALIGNMENT_KWARGS,
        ),
        sim,
        atol=tol,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
@pytest.mark.parametrize("skew", [False, True], ids=["normal", "skew"])
@pytest.mark.parametrize("m", [1, 2, 3, 4])
def test_ExactCFbend_multipoles(sim, m, skew):
    """Each curvilinear multipole of index m (quadrupole: m = 1), with the same
    field at r = 5 cm for all orders
    """
    k = [0.0] * (m + 1)
    k[m] = factorial(m) * 5.0e-3 / 0.05**m
    zero = [0.0] * (m + 1)
    el = elements.ExactCFbend(
        ds=0.15,
        k_normal=[2.0] + (zero if skew else k)[1:],
        k_skew=k if skew else zero,
    )

    tol = atol
    if Config.precision == "SINGLE":
        # FIXME: the bend map cancels terms of order rc in single precision,
        # M^T J M - J is up to 5e-4; with #1774 and #1776, it is up to 3e-5 and
        # the default atol holds
        # https://github.com/BLAST-ImpactX/impactx/issues/1771
        tol = 1.0e-3
    check_symplectic(el, sim, atol=tol)


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
@pytest.mark.parametrize(("unit", "k"), [(0, 1.0), (1, 3.5)], ids=["madx", "marylie"])
def test_ExactQuad(sim, unit, k):
    check_symplectic(
        elements.ExactQuad(
            ds=1.0,
            k=k,
            unit=unit,
            int_order=4,
            nslice=nslice,
            mapsteps=mapsteps,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


def test_ExactSbend(sim):
    tol = atol
    if Config.precision == "SINGLE":
        # FIXME: the bend map cancels terms of order rc in single precision,
        # M^T J M - J is up to 6e-4; with #1776, it is up to 4e-6 and the default
        # atol holds
        # https://github.com/BLAST-ImpactX/impactx/issues/1771
        tol = 1.2e-3
    check_symplectic(
        elements.ExactSbend(
            ds=1.0, phi=10.0, B=0.45, nslice=nslice, **ALIGNMENT_KWARGS
        ),
        sim,
        atol=tol,
    )


def test_Quad(sim):
    check_symplectic(
        elements.Quad(ds=1.0, k=1.0, nslice=nslice, **ALIGNMENT_KWARGS), sim
    )


def test_Sbend(sim):
    check_symplectic(
        elements.Sbend(ds=0.5, rc=-10.346, nslice=nslice, **ALIGNMENT_KWARGS), sim
    )


def test_SoftQuadrupole(sim):
    check_symplectic(
        elements.SoftQuadrupole(
            ds=1.0,
            gscale=1.0,
            cos_coefficients=[2.0, 0.3, -0.2, 0.05],
            sin_coefficients=[0.0, 0.25, -0.1, 0.04],
            mapsteps=mapsteps,
            nslice=nslice,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


@pytest.mark.parametrize(
    ("unit", "bscale"),
    [(0, 1.233482899483985), (1, 2.0)],
    ids=["madx", "si"],
)
def test_SoftSolenoid(sim, unit, bscale):
    check_symplectic(
        elements.SoftSolenoid(
            ds=6.0,
            bscale=bscale,
            cos_coefficients=[
                0.350807812299706,
                0.323554693720069,
                0.260320578919415,
                0.182848575294969,
                0.106921016050403,
                4.409581845710694e-002,
                -9.416427163897508e-006,
                -2.459452716865687e-002,
                -3.272762575737291e-002,
                -2.936414401076162e-002,
                -1.995780078926890e-002,
                -9.102893342953847e-003,
                -2.456410658713271e-006,
                5.788233017324325e-003,
                8.040408292420691e-003,
                7.480064552867431e-003,
                5.230254569468851e-003,
                2.447614547094685e-003,
                -1.095525090532255e-006,
                -1.614586867387170e-003,
                -2.281365457438345e-003,
                -2.148709081338292e-003,
                -1.522541739363011e-003,
                -7.185505862719508e-004,
                -6.171194824600157e-007,
                4.842109305036943e-004,
                6.874508102002901e-004,
                6.535550288205728e-004,
                4.648795813759210e-004,
                2.216564722797528e-004,
                -4.100982995210341e-007,
                -1.499332112463395e-004,
                -2.151538438342482e-004,
                -2.044590946652016e-004,
                -1.468242784844341e-004,
            ],
            sin_coefficients=[0.0, 0.1, -0.03, 0.02, -0.01, 0.006] + [0.0] * 29,
            unit=unit,
            mapsteps=mapsteps,
            nslice=nslice,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


def test_Sol(sim):
    check_symplectic(
        elements.Sol(ds=3.820395, ks=0.8223219329893234, **ALIGNMENT_KWARGS), sim
    )


# =============================================================================
# Thin kick elements
# =============================================================================


def test_Buncher(sim):
    check_symplectic(elements.Buncher(V=0.01, k=15.0, **ALIGNMENT_KWARGS), sim)


@pytest.mark.parametrize(
    ("unit", "xkick", "ykick"),
    [
        ("dimensionless", 2.0e-3, 3.0e-3),
        ("T-m", 7.0e-3, -4.0e-3),
    ],
    ids=["dimensionless", "tm"],
)
def test_Kicker(sim, unit, xkick, ykick):
    check_symplectic(
        elements.Kicker(xkick=xkick, ykick=ykick, unit=unit, **ALIGNMENT_KWARGS), sim
    )


@pytest.mark.parametrize(
    ("multipole", "K_normal", "K_skew"),
    [(4, 65.0, 6.0), (2, 300.0, 100.0)],
    ids=["octupole", "strong_quadrupole"],
)
def test_Multipole(sim, multipole, K_normal, K_skew):
    tol = atol
    if Config.precision == "SINGLE" and multipole == 2:
        # the strong kick makes the Jacobian entries of order 300, and the rounding
        # error of M^T J M grows with them, up to 7e-4
        tol = 1.5e-3
    check_symplectic(
        elements.Multipole(
            multipole=multipole, K_normal=K_normal, K_skew=K_skew, **ALIGNMENT_KWARGS
        ),
        sim,
        atol=tol,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_NonlinearLens(sim):
    check_symplectic(
        elements.NonlinearLens(knll=4.0e-6, cnll=0.01, **ALIGNMENT_KWARGS), sim
    )


@pytest.mark.parametrize(("unit", "k"), [(0, 1.0 / 0.5), (1, 6.0)], ids=["madx", "si"])
def test_TaperedPL(sim, unit, k):
    check_symplectic(
        elements.TaperedPL(
            k=k, taper=11.488289081903567, unit=unit, **ALIGNMENT_KWARGS
        ),
        sim,
    )


def test_ThinDipole(sim):
    check_symplectic(elements.ThinDipole(theta=0.45, rc=1.0, **ALIGNMENT_KWARGS), sim)


# =============================================================================
# Rotation elements
# =============================================================================


def test_PlaneXYRot(sim):
    check_symplectic(elements.PlaneXYRot(angle=90.0, **ALIGNMENT_KWARGS), sim)


def test_PRot(sim):
    check_symplectic(elements.PRot(phi_in=0.0, phi_out=-5.0), sim)


# =============================================================================
# RF / energy-changing elements
# =============================================================================


def test_ShortRF(sim):
    check_symplectic(
        elements.ShortRF(V=1000.0, freq=1.3e9, phase=-89.5, **ALIGNMENT_KWARGS), sim
    )


def test_RFCavity(sim):
    check_symplectic(
        elements.RFCavity(
            ds=1.31879807,
            escale=62.0,
            freq=1.3e9,
            phase=85.5,
            cos_coefficients=[
                0.1644024074311037,
                -0.1324009958969339,
                4.3443060026047219e-002,
                8.5602654094946495e-002,
                -0.2433578169042885,
                0.5297150596779437,
                0.7164884680963959,
                -5.2579522442877296e-003,
                -5.5025369142193678e-002,
                4.6845673335028933e-002,
                -2.3279346335638568e-002,
                4.0800777539657775e-003,
                4.1378326533752169e-003,
                -2.5040533340490805e-003,
                -4.0654981400000964e-003,
                9.6630592067498289e-003,
                -8.5275895985990214e-003,
                -5.8078747006425020e-002,
                -2.4044337836660403e-002,
                1.0968240064697212e-002,
                -3.4461179858301418e-003,
                -8.1201564869443749e-004,
                2.1438992904959380e-003,
                -1.4997753525697276e-003,
                1.8685171825676386e-004,
            ],
            sin_coefficients=[0.0, 0.1, -0.03, 0.02, -0.01] + [0.0] * 20,
            mapsteps=mapsteps,
            nslice=nslice,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


# =============================================================================
# Edge elements
# =============================================================================


@pytest.mark.parametrize("location", ["entry", "exit"])
@pytest.mark.parametrize(
    ("model", "g", "modify_ref_part"),
    [
        ("linear_slac75", 0.0, False),
        ("linear_slac75", 0.058, False),
        ("linear_hwang_lee", 0.058, False),
        ("nonlinear_hwang_lee", 0.058, False),
        ("nonlinear_hwang_lee", 0.058, True),
    ],
    ids=[
        "slac75-zero-gap",
        "slac75",
        "hwang-lee-linear",
        "hwang-lee-nonlinear",
        "hwang-lee-nonlinear-modify-ref",
    ],
)
def test_DipEdge(sim, model, g, modify_ref_part, location):
    check_symplectic(
        elements.DipEdge(
            psi=-0.048345620280243,
            rc=-10.3462283686195526,
            g=g,
            # nonzero fringe field integrals, for all terms of the Hwang-Lee map
            K1=0.3,
            K2=0.5,
            K4=0.2,
            K5=0.1,
            K6=0.05,
            model=model,
            location=location,
            modify_ref_part=modify_ref_part,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


@pytest.mark.parametrize(("unit", "k"), [(0, 1.0), (1, 3.5)], ids=["madx", "marylie"])
@pytest.mark.parametrize("flag", ["entry", "exit"])
def test_QuadEdge(sim, flag, unit, k):
    check_symplectic(
        elements.QuadEdge(k=k, unit=unit, flag=flag, **ALIGNMENT_KWARGS), sim
    )


# =============================================================================
# Linear map
# =============================================================================


def test_LinearMap(sim):
    R1 = Map6x6.identity()
    ds1 = 0.25
    ct = np.cos(0.31)
    st = np.sin(0.31)
    lt = ds1 / 16.6464  # ds / (beta*gamma^2)

    R1[1, 1] = ct
    R1[1, 2] = ct * ds1
    R1[1, 3] = st
    R1[1, 4] = st * ds1
    R1[2, 1] = 0.0
    R1[2, 2] = ct
    R1[2, 3] = 0.0
    R1[2, 4] = st
    R1[3, 1] = -st
    R1[3, 2] = -st * ds1
    R1[3, 3] = ct
    R1[3, 4] = ct * ds1
    R1[4, 1] = 0.0
    R1[4, 2] = -st
    R1[4, 3] = 0.0
    R1[4, 4] = ct
    R1[5, 6] = lt

    check_symplectic(elements.LinearMap(R=R1, ds=ds1, **ALIGNMENT_KWARGS), sim)
