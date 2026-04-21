#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# This set of tests verifies that forward + inverse stepping through
# each element returns particles (and spin) to their initial state.
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import (
    Config,
    ImpactX,
    Map6x6,
    distribution,
    elements,
    twiss,
)

# test config
npart = 1_000

# element config
nslice = 1
mapsteps = 4

ALIGNMENT_KWARGS = dict(dx=1.0e-3, dy=-2.0e-3, rotation=0.5)
PIPE_KWARGS = dict(ALIGNMENT_KWARGS, aperture_x=0.03, aperture_y=0.04)

if Config.precision == "SINGLE":
    phase_atol = 2.0e-6
    spin_atol = 3.0e-7
    ref_atol = 3.0e-4
else:
    phase_atol = 1.0e-10
    spin_atol = 1.0e-12
    ref_atol = 1.0e-10


@pytest.fixture(scope="function")
def sim(request):
    class SimContextManager:
        def __enter__(self):
            self.sim = ImpactX()

            # set numerical parameters and IO control
            self.sim.particle_shape = 2  # B-spline order
            self.sim.space_charge = False
            self.sim.diagnostics = False
            self.sim.slice_step_diagnostics = False

            spin = getattr(
                request, "param", False
            )  # default to False if not parametrized
            self.sim.spin = spin

            self.sim.init_grids()

            beam = self.sim.beam
            beam.clear_particles()

            # load a 1 GeV electron beam with an initial
            # unnormalized rms emittance of 1 nm
            kin_energy_MeV = 1.0e3  # reference energy
            bunch_charge_C = 1.0e-9  # used with space charge

            beam.ref.set_species("electron").set_kin_energy_MeV(kin_energy_MeV)

            #   particle bunch
            distr = distribution.Waterbag(
                **twiss(
                    beta_x=1.0,
                    beta_y=1.0,
                    beta_t=1.0,
                    emitt_x=1e-09,
                    emitt_y=1e-09,
                    emitt_t=1e-06,
                    alpha_x=0.0,
                    alpha_y=0.0,
                    alpha_t=0.0,
                )
            )
            if spin:
                spin_vectors = distribution.SpinvMF(
                    0.6,
                    0.5,
                    0.4,
                )
                self.sim.add_particles(bunch_charge_C, distr, npart, spin_vectors)
            else:
                self.sim.add_particles(bunch_charge_C, distr, npart)
            assert self.sim.beam.total_number_of_particles() == npart

            return self.sim

        def __exit__(self, exc_type, exc_value, traceback):
            self.sim.finalize()
            del self.sim

    with SimContextManager() as sim:
        yield sim


PHASE_COLS = [
    "position_x",
    "position_y",
    "position_t",
    "momentum_x",
    "momentum_y",
    "momentum_t",
]
SPIN_COLS = ["spin_x", "spin_y", "spin_z"]
REF_COLS = ["x", "y", "z", "t", "px", "py", "pz", "pt", "s"]


def save_state(beam, spin=False):
    """Save a copy of particle phase space (and optionally spin) arrays."""
    df = beam.to_df()
    cols = PHASE_COLS + (SPIN_COLS if spin else [])
    return {c: df[c].to_numpy().copy() for c in cols}


def save_ref_state(beam):
    """Save a copy of the reference particle state."""
    return beam.ref.copy()


def check_changed(beam, initial, spin=False):
    """Assert that at least one coordinate changed meaningfully."""
    df = beam.to_df()
    cols = PHASE_COLS + (SPIN_COLS if spin else [])
    any_changed = False
    for c in cols:
        if not np.allclose(df[c].to_numpy(), initial[c], atol=1e-14, rtol=0):
            any_changed = True
            break
    assert any_changed, "Forward push did not change any particle coordinate"


def check_roundtrip(
    beam,
    initial,
    initial_ref,
    phase_atol=phase_atol,
    spin_atol=spin_atol,
    ref_atol=ref_atol,
    spin=False,
):
    """Assert particles and the reference particle returned to initial state."""
    df = beam.to_df()
    for c in PHASE_COLS:
        np.testing.assert_allclose(
            df[c].to_numpy(),
            initial[c],
            atol=phase_atol,
            rtol=0,
            err_msg=f"Roundtrip mismatch in {c}",
        )
    if spin:
        for c in SPIN_COLS:
            np.testing.assert_allclose(
                df[c].to_numpy(),
                initial[c],
                atol=spin_atol,
                rtol=0,
                err_msg=f"Roundtrip mismatch in {c}",
            )
    for c in REF_COLS:
        np.testing.assert_allclose(
            getattr(beam.ref, c),
            getattr(initial_ref, c),
            atol=ref_atol,
            rtol=0,
            err_msg=f"Reference-particle roundtrip mismatch in {c}",
        )


def roundtrip(
    el,
    sim,
    phase_atol=phase_atol,
    spin_atol=spin_atol,
    ref_atol=ref_atol,
    spin=False,
):
    """Run forward + reverse + forward and verify roundtrip."""
    beam = sim.beam

    initial = save_state(beam, spin=spin)
    initial_ref = save_ref_state(beam)

    # forward push
    el.push(beam)

    # verify something changed
    check_changed(beam, initial, spin=spin)

    # reverse and push again
    el.reverse()
    el.push(beam)

    # verify roundtrip
    check_roundtrip(
        beam,
        initial,
        initial_ref=initial_ref,
        phase_atol=phase_atol,
        spin_atol=spin_atol,
        ref_atol=ref_atol,
        spin=spin,
    )


# =============================================================================
# Thick elements (negate ds)
# =============================================================================


# CFbend is spin-capable (inherits mixin::SpinTransport) but its spin
# reversibility path does not round-trip at machine precision
# (https://github.com/BLAST-ImpactX/impactx/issues/1328). Keep nospin-only
# until the C++ spin-reversal is fixed.
def test_CFbend(sim):
    roundtrip(
        elements.CFbend(
            ds=0.5,
            rc=7.613657587094493,
            k=-7.057403,
            nslice=nslice,
            **PIPE_KWARGS,
        ),
        sim,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_ChrDrift(sim):
    roundtrip(
        elements.ChrDrift(ds=0.25, nslice=nslice, **PIPE_KWARGS),
        sim,
        spin=sim.spin,
    )


@pytest.mark.parametrize(
    ("unit", "k"),
    [(0, 2.98636067687944129), (1, 5.0)],
    ids=["madx", "si"],
)
def test_ChrPlasmaLens(sim, unit, k):
    roundtrip(
        elements.ChrPlasmaLens(
            ds=0.331817852986604588,
            k=k,
            unit=unit,
            nslice=nslice,
            **PIPE_KWARGS,
        ),
        sim,
        phase_atol=2e-6 if unit == 1 else phase_atol,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_ChrQuad(sim):
    kwargs = {} if sim.spin else PIPE_KWARGS
    roundtrip(
        elements.ChrQuad(ds=1.0, k=1.0, unit=0, nslice=nslice, **kwargs),
        sim,
        spin=sim.spin,
    )


def test_ChrQuad_marylie_unit(sim):
    roundtrip(
        elements.ChrQuad(ds=1.0, k=3.5, unit=1, nslice=nslice, **PIPE_KWARGS),
        sim,
    )


@pytest.mark.parametrize("sim", [True], indirect=True, ids=["spin"])
def test_ChrQuad_marylie_unit_spin(sim):
    roundtrip(
        elements.ChrQuad(ds=1.0, k=1.0, unit=1, nslice=nslice),
        sim,
        spin=sim.spin,
    )


def test_ChrAcc(sim):
    roundtrip(
        elements.ChrAcc(
            ds=1.8,
            ez=10871.950994502130424,
            bz=250.0,
            nslice=nslice,
            **PIPE_KWARGS,
        ),
        sim,
    )


def test_ConstF(sim):
    roundtrip(
        elements.ConstF(ds=2.0, kx=1.0, ky=1.0, kt=1.0, nslice=nslice, **PIPE_KWARGS),
        sim,
        # ConstF's longitudinal (t, pt) block is the SP-sensitive outlier in this file.
        phase_atol=1.0e-3 if Config.precision == "SINGLE" else phase_atol,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_Drift(sim):
    roundtrip(
        elements.Drift(ds=0.25, nslice=nslice, **PIPE_KWARGS),
        sim,
        spin=sim.spin,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_ExactDrift(sim):
    roundtrip(
        elements.ExactDrift(ds=0.25, nslice=nslice, **PIPE_KWARGS),
        sim,
        spin=sim.spin,
    )


@pytest.mark.parametrize(
    ("unit", "k_normal", "k_skew"),
    [
        (0, [0.0, 1.0, 0.15], [0.0, 0.35, -0.05]),
        (1, [0.0, 3.5, 0.45], [0.0, 0.8, -0.12]),
    ],
    ids=["madx", "si"],
)
def test_ExactMultipole(sim, unit, k_normal, k_skew):
    roundtrip(
        elements.ExactMultipole(
            ds=1.0,
            k_normal=k_normal,
            k_skew=k_skew,
            unit=unit,
            mapsteps=mapsteps,
            nslice=nslice,
            int_order=4,
            **PIPE_KWARGS,
        ),
        sim,
        phase_atol=1e-8,
    )


@pytest.mark.parametrize(
    ("unit", "k_normal", "k_skew"),
    [
        (0, [0.12, 0.015, -0.002], [0.0, 0.01, 0.001]),
        (1, [0.4, 0.05, -0.008], [0.0, 0.03, 0.004]),
    ],
    ids=["madx", "si"],
)
def test_ExactCFbend(sim, unit, k_normal, k_skew):
    roundtrip(
        elements.ExactCFbend(
            ds=1.0,
            k_normal=k_normal,
            k_skew=k_skew,
            unit=unit,
            int_order=4,
            mapsteps=mapsteps,
            nslice=nslice,
            **PIPE_KWARGS,
        ),
        sim,
        phase_atol=1e-8,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
@pytest.mark.parametrize(("unit", "k"), [(0, 1.0), (1, 3.5)], ids=["madx", "marylie"])
def test_ExactQuad(sim, unit, k):
    kwargs = {} if sim.spin else PIPE_KWARGS
    roundtrip(
        elements.ExactQuad(
            ds=1.0,
            k=k,
            unit=unit,
            int_order=4,
            nslice=nslice,
            mapsteps=mapsteps,
            **kwargs,
        ),
        sim,
        phase_atol=1e-8,
        spin=sim.spin,
    )


def test_ExactSbend(sim):
    roundtrip(
        elements.ExactSbend(ds=1.0, phi=10.0, B=0.45, nslice=nslice, **PIPE_KWARGS),
        sim,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_Quad(sim):
    kwargs = {} if sim.spin else PIPE_KWARGS
    roundtrip(
        elements.Quad(ds=1.0, k=1.0, nslice=nslice, **kwargs),
        sim,
        spin=sim.spin,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_Sbend(sim):
    kwargs = {} if sim.spin else PIPE_KWARGS
    roundtrip(
        elements.Sbend(ds=0.5, rc=-10.346, nslice=nslice, **kwargs),
        sim,
        spin=sim.spin,
    )


def test_SoftQuadrupole(sim):
    roundtrip(
        elements.SoftQuadrupole(
            ds=1.0,
            gscale=1.0,
            cos_coefficients=[2.0, 0.3, -0.2, 0.05],
            sin_coefficients=[0.0, 0.25, -0.1, 0.04],
            mapsteps=mapsteps,
            nslice=nslice,
            **PIPE_KWARGS,
        ),
        sim,
        phase_atol=1e-8,
    )


@pytest.mark.parametrize(
    ("unit", "bscale"),
    [(0, 1.233482899483985), (1, 2.0)],
    ids=["madx", "si"],
)
def test_SoftSolenoid(sim, unit, bscale):
    roundtrip(
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
            **PIPE_KWARGS,
        ),
        sim,
        phase_atol=1e-8,
    )


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
def test_Sol(sim):
    kwargs = {} if sim.spin else PIPE_KWARGS
    roundtrip(
        elements.Sol(ds=3.820395, ks=0.8223219329893234, **kwargs),
        sim,
        spin=sim.spin,
    )


# =============================================================================
# Thin kick elements (negate strength)
# =============================================================================


def test_Buncher(sim):
    roundtrip(elements.Buncher(V=0.01, k=15.0, **ALIGNMENT_KWARGS), sim)


@pytest.mark.parametrize(
    ("unit", "xkick", "ykick"),
    [
        ("dimensionless", 2.0e-3, 3.0e-3),
        ("T-m", 7.0e-3, -4.0e-3),
    ],
    ids=["dimensionless", "tm"],
)
def test_Kicker(sim, unit, xkick, ykick):
    roundtrip(
        elements.Kicker(xkick=xkick, ykick=ykick, unit=unit, **ALIGNMENT_KWARGS),
        sim,
    )


def test_Multipole(sim):
    roundtrip(
        elements.Multipole(multipole=4, K_normal=65.0, K_skew=6.0, **ALIGNMENT_KWARGS),
        sim,
    )


def test_NonlinearLens(sim):
    roundtrip(elements.NonlinearLens(knll=4.0e-6, cnll=0.01, **ALIGNMENT_KWARGS), sim)


@pytest.mark.parametrize(("unit", "k"), [(0, 1.0 / 0.5), (1, 6.0)], ids=["madx", "si"])
def test_TaperedPL(sim, unit, k):
    roundtrip(
        elements.TaperedPL(
            k=k,
            taper=11.488289081903567,
            unit=unit,
            **ALIGNMENT_KWARGS,
        ),
        sim,
    )


def test_ThinDipole(sim):
    roundtrip(elements.ThinDipole(theta=0.45, rc=1.0, **ALIGNMENT_KWARGS), sim)


# =============================================================================
# Rotation elements
# =============================================================================


def test_PlaneXYRot(sim):
    roundtrip(elements.PlaneXYRot(angle=90.0, **ALIGNMENT_KWARGS), sim)


def test_PRot(sim):
    roundtrip(elements.PRot(phi_in=0.0, phi_out=-5.0), sim)


# =============================================================================
# RF / energy-changing elements
# =============================================================================


def test_ShortRF(sim):
    roundtrip(
        elements.ShortRF(V=1000.0, freq=1.3e9, phase=-89.5, **ALIGNMENT_KWARGS), sim
    )


def test_RFCavity(sim):
    roundtrip(
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
            **PIPE_KWARGS,
        ),
        sim,
        phase_atol=5e-4,  # energy-changing element: finite mapsteps limit roundtrip precision
    )


# =============================================================================
# Edge elements
# =============================================================================


@pytest.mark.parametrize(
    ("model", "g", "K2"),
    [
        ("linear", 0.0, 0.0),
        ("linear", 0.058, 0.5),
        ("nonlinear", 0.058, 0.5),
    ],
    ids=["linear-zero-gap", "linear-finite-gap", "nonlinear-finite-gap"],
)
def test_DipEdge(sim, model, g, K2):
    rc = 10.3462283686195526
    psi = 0.048345620280243
    roundtrip(
        elements.DipEdge(
            psi=-psi,
            rc=-rc,
            g=g,
            K2=K2,
            model=model,
            **ALIGNMENT_KWARGS,
        ),
        sim,
        phase_atol=2e-9 if model == "nonlinear" else phase_atol,
    )


@pytest.mark.parametrize(("unit", "k"), [(0, 1.0), (1, 3.5)], ids=["madx", "marylie"])
@pytest.mark.parametrize("flag", ["entry", "exit"])
def test_QuadEdge(sim, flag, unit, k):
    roundtrip(
        elements.QuadEdge(k=k, unit=unit, flag=flag, **ALIGNMENT_KWARGS),
        sim,
    )


# =============================================================================
# Linear map
# =============================================================================


@pytest.mark.parametrize("sim", [True, False], indirect=True, ids=["spin", "nospin"])
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

    el = elements.LinearMap(R=R1, ds=ds1, **ALIGNMENT_KWARGS)
    assert el.symplectic
    roundtrip(el, sim, spin=sim.spin)


# There is a hacky implementation draft in the C++ code for
# reversibility of the SpinMap, but it needs a rework.
#
# @pytest.mark.parametrize("sim", [True], indirect=True, ids=["spin"])
# def test_SpinMap(sim):
#     v = Vector3()
#     A = Map3x6()
#
#     v[1] = 0.17
#     v[2] = -0.11
#     v[3] = 0.07
#
#     A[1, 1] = 0.21
#     A[1, 4] = -0.13
#     A[2, 2] = 0.19
#     A[2, 5] = 0.09
#     A[3, 3] = -0.15
#     A[3, 6] = 0.05
#
#    roundtrip(elements.SpinMap(v=v, A=A), sim, spin=sim.spin)
