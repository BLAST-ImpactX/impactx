from __future__ import annotations

import typing

import amrex.space3d.amrex_3d_pybind
import impactx.impactx_pybind

__all__: list[str] = [
    "alpha",
    "convolve_fft",
    "deposit_charge",
    "derivative_charge",
    "unit_step",
    "w_l_csr",
    "w_l_rf",
    "w_t_rf",
]

def alpha(arg0: typing.SupportsFloat) -> float:
    """
    Alpha Function
    """

def convolve_fft(
    arg0: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg1: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg2: typing.SupportsFloat,
) -> amrex.space3d.amrex_3d_pybind.PODVector_real_std:
    """
    FFT Convolution
    """

def deposit_charge(
    arg0: impactx.impactx_pybind.ImpactXParticleContainer,
    arg1: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg2: typing.SupportsFloat,
    arg3: typing.SupportsFloat,
    arg4: bool,
) -> None:
    """
    Deposit Charge Distribution Function
    """

def derivative_charge(
    arg0: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg1: amrex.space3d.amrex_3d_pybind.PODVector_real_std,
    arg2: typing.SupportsFloat,
    arg3: bool,
) -> None:
    """
    Derivative of Charge Profile Function
    """

def unit_step(arg0: typing.SupportsFloat) -> float:
    """
    Step Function
    """

def w_l_csr(
    arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat
) -> float:
    """
    CSR Wake Function
    """

def w_l_rf(
    arg0: typing.SupportsFloat,
    arg1: typing.SupportsFloat,
    arg2: typing.SupportsFloat,
    arg3: typing.SupportsFloat,
) -> float:
    """
    Longitudinal Resistive Wall Wake Function
    """

def w_t_rf(
    arg0: typing.SupportsFloat,
    arg1: typing.SupportsFloat,
    arg2: typing.SupportsFloat,
    arg3: typing.SupportsFloat,
) -> float:
    """
    Transverse Resistive Wall Wake Function
    """
