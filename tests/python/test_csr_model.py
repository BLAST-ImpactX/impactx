#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math

import numpy as np
import pytest
from scipy.constants import c

from impactx import Config, ImpactX, distribution, elements

# beam parameters, following examples/chicane/run_chicane_csr.py
KIN_ENERGY_MEV = 5.0e3
BUNCH_CHARGE_C = 1.0e-9
NPART = 2000

BUILD_DTYPE = np.float64 if Config.precision == "DOUBLE" else np.float32
CONTEXT_RTOL = 1.0e-14 if Config.precision == "DOUBLE" else 2.0e-6
SLICE_POSITION_RTOL = 1.0e-13 if Config.precision == "DOUBLE" else 2.0e-6
REFERENCE_RTOL = 1.0e-12 if Config.precision == "DOUBLE" else 2.0e-6
CHARGE_RTOL = 1.0e-6 if Config.precision == "DOUBLE" else 2.0e-5
KICK_RTOL = 1.0e-10 if Config.precision == "DOUBLE" else 2.0e-3
ZERO_ATOL = 1.0e-14


def force_for_precision():
    """A force large enough to resolve its momentum change at build precision."""
    return 1.0e-15 if Config.precision == "DOUBLE" else 1.0e-12


def build_sim(
    csr=True,
    csr_bins=32,
    ds=0.5,
    rc=10.35,
    nslice=6,
    kick_model=None,
    distr=None,
):
    """Build a one-bend simulation for CSR kick model tests."""
    sim = ImpactX()

    sim.particle_shape = 2
    sim.space_charge = False
    sim.csr = csr
    sim.csr_bins = csr_bins
    sim.diagnostics = False
    sim.slice_step_diagnostics = False
    if kick_model is not None:
        sim.csr_kick_model = kick_model

    sim.init_grids()

    ref = sim.beam.ref
    ref.set_species("electron").set_kin_energy_MeV(KIN_ENERGY_MEV)

    if distr is None:
        distr = distribution.Gaussian(
            lambdaX=2.2951017632e-5,
            lambdaY=1.3084093142e-5,
            lambdaT=5.5555553e-8,
            lambdaPx=1.598353425e-6,
            lambdaPy=2.803697378e-6,
            lambdaPt=2.000000000e-6,
            muxpx=0.933345606203060,
            muypy=0.933345606203060,
            mutpt=0.999999961419755,
        )
    sim.add_particles(BUNCH_CHARGE_C, distr, NPART)

    sim.lattice.append(elements.Sbend(name="bend1", ds=ds, rc=rc, nslice=nslice))

    return sim


def ref_momentum_SI(ref):
    """Reference momentum p_ref = beta * gamma * m * c in SI [kg m/s]."""
    return ref.beta_gamma * ref.mass * c


def mean_momenta(sim):
    """Mean momenta (px, py, pt) of the beam (dimensionless, p/p_ref)."""
    df = sim.beam.to_df()
    return (
        df["momentum_x"].mean(),
        df["momentum_y"].mean(),
        df["momentum_t"].mean(),
    )


def test_csr_kick_model_context_and_profile():
    """A zero-kick model receives the correct context and profile and does
    not alter the beam momentum."""
    ds = 0.5
    rc = -10.35
    nslice = 6
    csr_bins = 32
    calls = []

    def recorder(profile, context):
        charge = profile.charge.to_numpy(copy=True)
        mean_x = profile.mean_x.to_numpy(copy=True)
        mean_y = profile.mean_y.to_numpy(copy=True)
        calls.append(
            {
                "element_name": context.element_name,
                "element_type": context.element_type,
                "rc": context.rc,
                "signed_rc": context.signed_rc,
                "ds": context.ds,
                "nslice": context.nslice,
                "slice": context.slice,
                "s": context.s,
                "slice_ds": context.slice_ds,
                "beta_gamma": context.ref.beta_gamma,
                "num_bins": profile.num_bins,
                "bin_min": profile.bin_min,
                "bin_size": profile.bin_size,
                "charge": charge,
                "mean_x": mean_x,
                "mean_y": mean_y,
            }
        )
        return np.zeros(profile.num_bins)

    sim = build_sim(csr_bins=csr_bins, ds=ds, rc=rc, nslice=nslice, kick_model=recorder)
    try:
        beta_gamma_ref = sim.beam.ref.beta_gamma
        _, _, pt_mean_initial = mean_momenta(sim)

        sim.track_particles()

        # called once per slice, in order
        assert len(calls) == nslice
        assert [call["slice"] for call in calls] == list(range(nslice))

        for call in calls:
            # element context
            assert call["element_name"] == "bend1"
            assert call["element_type"] == "Sbend"
            assert math.isclose(call["rc"], abs(rc), rel_tol=CONTEXT_RTOL)
            assert math.isclose(call["signed_rc"], rc, rel_tol=CONTEXT_RTOL)
            assert math.isclose(call["ds"], ds, rel_tol=CONTEXT_RTOL)
            assert call["nslice"] == nslice
            assert math.isclose(call["slice_ds"], ds / nslice, rel_tol=CONTEXT_RTOL)
            assert math.isclose(
                call["s"],
                call["slice"] * call["slice_ds"],
                rel_tol=SLICE_POSITION_RTOL,
                abs_tol=0.0,
            ) or (call["slice"] == 0 and call["s"] == 0.0)
            assert math.isclose(
                call["beta_gamma"], beta_gamma_ref, rel_tol=REFERENCE_RTOL
            )

            # binned profile
            assert call["num_bins"] == csr_bins
            assert len(call["charge"]) == csr_bins + 1
            assert len(call["mean_x"]) == csr_bins + 1
            assert len(call["mean_y"]) == csr_bins + 1
            assert call["bin_size"] > 0.0
            # the histogram integrates to the bunch charge
            assert math.isclose(
                np.sum(call["charge"]) * call["bin_size"],
                BUNCH_CHARGE_C,
                rel_tol=CHARGE_RTOL,
            )
            assert np.all(np.isfinite(call["mean_x"]))
            assert np.all(np.isfinite(call["mean_y"]))

        # a zero kick does not alter the mean longitudinal momentum
        # (the bend map conserves pt)
        _, _, pt_mean_final = mean_momenta(sim)
        assert math.isclose(
            pt_mean_final, pt_mean_initial, rel_tol=0.0, abs_tol=ZERO_ATOL
        )
    finally:
        sim.finalize()


def run_constant_kick(force, as_dict):
    """Track through one bend with a constant longitudinal kick force [N]
    and return (mean pt change, ds, reference momentum [kg m/s])."""
    ds = 0.5

    def constant_kick(profile, context):
        kick = np.full(profile.num_bins, force)
        return {"pt": kick} if as_dict else kick

    sim = build_sim(ds=ds, kick_model=constant_kick)
    try:
        p_ref = ref_momentum_SI(sim.beam.ref)
        _, _, pt_mean_initial = mean_momenta(sim)
        sim.track_particles()
        _, _, pt_mean_final = mean_momenta(sim)
    finally:
        sim.finalize()

    return pt_mean_final - pt_mean_initial, ds, p_ref


def test_csr_kick_model_longitudinal_kick():
    """A constant longitudinal force F over a bend of length ds changes the
    beam mean pt by exactly -ds * F / (c * p_ref)."""
    force = force_for_precision()  # [N]
    delta_pt, ds, p_ref = run_constant_kick(force, as_dict=True)

    delta_pt_expected = -ds * force / (c * p_ref)
    assert math.isclose(delta_pt, delta_pt_expected, rel_tol=KICK_RTOL)


def test_csr_kick_model_bare_array_return():
    """Returning a bare array is equivalent to returning {'pt': array}."""
    force = force_for_precision()  # [N]
    delta_pt_dict, _, _ = run_constant_kick(force, as_dict=True)
    delta_pt_bare, _, _ = run_constant_kick(force, as_dict=False)

    assert math.isclose(
        delta_pt_dict,
        delta_pt_bare,
        rel_tol=KICK_RTOL,
        abs_tol=0.0,
    )


def test_csr_kick_model_transverse_kick():
    """A constant horizontal force F over a thin bend changes the beam mean
    px by ds * F / (beta * c * p_ref)."""
    ds = 1.0e-3
    force = 1.0e-12 if Config.precision == "DOUBLE" else 1.0e-10  # [N]

    def transverse_kick(profile, context):
        return {
            "pt": np.zeros(profile.num_bins),
            "px": np.full(profile.num_bins, force),
        }

    # an unchirped beam with negligible momentum spread, so that the bend
    # map's dispersion (px += -sin(theta) * pt) does not shift the mean px
    # (the default chicane beam is strongly chirped: sample <pt> ~ 5e-4)
    distr = distribution.Gaussian(
        lambdaX=2.2951017632e-5,
        lambdaY=1.3084093142e-5,
        lambdaT=5.5555553e-8,
        lambdaPx=1.598353425e-6,
        lambdaPy=2.803697378e-6,
        lambdaPt=1.0e-9,
        muxpx=0.0,
        muypy=0.0,
        mutpt=0.0,
    )
    sim = build_sim(ds=ds, rc=10.0, nslice=1, kick_model=transverse_kick, distr=distr)
    try:
        ref = sim.beam.ref
        p_ref = ref_momentum_SI(ref)
        beta_ref = ref.beta
        px_mean_initial, py_mean_initial, _ = mean_momenta(sim)
        sim.track_particles()
        px_mean_final, py_mean_final, _ = mean_momenta(sim)
    finally:
        sim.finalize()

    delta_px_expected = ds * force / (beta_ref * c * p_ref)
    # the bend map after the kick mixes phase space slightly (arc angle 1e-4)
    assert math.isclose(
        px_mean_final - px_mean_initial,
        delta_px_expected,
        rel_tol=1.0e-4 if Config.precision == "DOUBLE" else 2.0e-3,
    )
    # no vertical kick was applied
    assert math.isclose(py_mean_final, py_mean_initial, rel_tol=0.0, abs_tol=ZERO_ATOL)


def test_csr_kick_model_gpu_cupy_return():
    """On GPU builds, cupy arrays work as zero-copy inputs and as returns
    (via the CUDA array interface), including producer stream synchronization."""
    if Config.gpu_backend not in {"CUDA", "HIP"}:
        pytest.skip("CUDA/HIP-only test")
    cp = pytest.importorskip("cupy")

    ds = 0.5
    force = force_for_precision()  # [N]

    class ArrayWithStream:
        """CUDA array-interface wrapper that records the producer stream."""

        def __init__(self, array, stream):
            self.array = array
            self.stream = stream
            self.__cuda_array_interface__ = dict(array.__cuda_array_interface__)
            self.__cuda_array_interface__["stream"] = stream.ptr

    def cupy_kick(profile, context):
        lam = profile.charge.to_xp()  # zero-copy device view
        assert isinstance(lam, cp.ndarray)
        assert lam.shape == (profile.num_bins + 1,)
        stream = cp.cuda.Stream(non_blocking=True)
        with stream:
            kick = cp.full(profile.num_bins, force, dtype=BUILD_DTYPE)
        return {"pt": ArrayWithStream(kick, stream)}

    sim = build_sim(ds=ds, kick_model=cupy_kick)
    try:
        p_ref = ref_momentum_SI(sim.beam.ref)
        _, _, pt_mean_initial = mean_momenta(sim)
        sim.track_particles()
        _, _, pt_mean_final = mean_momenta(sim)
    finally:
        sim.finalize()

    delta_pt_expected = -ds * force / (c * p_ref)
    assert math.isclose(
        pt_mean_final - pt_mean_initial, delta_pt_expected, rel_tol=KICK_RTOL
    )


def test_csr_kick_model_disabled_csr():
    """With csr disabled, a set kick model is never called."""
    calls = []

    def recorder(profile, context):
        calls.append(context.slice)
        return np.zeros(profile.num_bins)

    sim = build_sim(csr=False, kick_model=recorder)
    try:
        sim.track_particles()
    finally:
        sim.finalize()

    assert calls == []


@pytest.mark.parametrize(
    "bad_return, match",
    [
        (lambda n: np.zeros(n - 1), "wrong length"),
        (lambda n: {"px": np.zeros(n)}, "missing the required key 'pt'"),
        (lambda n: {"pt": np.zeros(n), "foo": np.zeros(n)}, "unknown key"),
        (lambda n: object(), "cannot be converted"),
    ],
)
def test_csr_kick_model_invalid_return(bad_return, match):
    """Invalid model return values raise errors naming the element."""

    def bad_model(profile, context):
        return bad_return(profile.num_bins)

    sim = build_sim(kick_model=bad_model)
    try:
        with pytest.raises(RuntimeError, match=match):
            sim.track_particles()
    finally:
        sim.finalize()


def test_csr_kick_model_property_validation():
    """The csr_kick_model property accepts callables and None only."""
    sim = ImpactX()
    try:
        assert sim.csr_kick_model is None

        def model(profile, context):
            return np.zeros(profile.num_bins)

        sim.csr_kick_model = model
        assert sim.csr_kick_model is model

        sim.csr_kick_model = None
        assert sim.csr_kick_model is None

        with pytest.raises(TypeError, match="callable"):
            sim.csr_kick_model = 5
    finally:
        sim.finalize()
