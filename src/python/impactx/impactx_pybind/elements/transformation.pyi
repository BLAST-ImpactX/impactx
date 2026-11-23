"""
Transform and modify lattices
"""

from __future__ import annotations

import typing

import impactx.impactx_pybind.elements

__all__: list[str] = ["insert_element_every_ds"]

def insert_element_every_ds(
    list: typing.Any,
    ds: typing.SupportsFloat | typing.SupportsIndex,
    element: impactx.impactx_pybind.elements.Empty
    | impactx.impactx_pybind.elements.Aperture
    | impactx.impactx_pybind.elements.Buncher
    | impactx.impactx_pybind.elements.CFbend
    | impactx.impactx_pybind.elements.ChrAcc
    | impactx.impactx_pybind.elements.ChrDrift
    | impactx.impactx_pybind.elements.ChrPlasmaLens
    | impactx.impactx_pybind.elements.ChrQuad
    | impactx.impactx_pybind.elements.ConstF
    | impactx.impactx_pybind.elements.BeamMonitor
    | impactx.impactx_pybind.elements.DipEdge
    | impactx.impactx_pybind.elements.Drift
    | impactx.impactx_pybind.elements.ExactCFbend
    | impactx.impactx_pybind.elements.ExactDrift
    | impactx.impactx_pybind.elements.ExactMultipole
    | impactx.impactx_pybind.elements.ExactQuad
    | impactx.impactx_pybind.elements.ExactSbend
    | impactx.impactx_pybind.elements.Kicker
    | impactx.impactx_pybind.elements.LinearMap
    | impactx.impactx_pybind.elements.Marker
    | impactx.impactx_pybind.elements.Multipole
    | impactx.impactx_pybind.elements.NonlinearLens
    | impactx.impactx_pybind.elements.PlaneXYRot
    | impactx.impactx_pybind.elements.PolygonAperture
    | impactx.impactx_pybind.elements.Programmable
    | impactx.impactx_pybind.elements.PRot
    | impactx.impactx_pybind.elements.Quad
    | impactx.impactx_pybind.elements.QuadEdge
    | impactx.impactx_pybind.elements.RFCavity
    | impactx.impactx_pybind.elements.Sbend
    | impactx.impactx_pybind.elements.ShortRF
    | impactx.impactx_pybind.elements.SoftSolenoid
    | impactx.impactx_pybind.elements.SoftQuadrupole
    | impactx.impactx_pybind.elements.Sol
    | impactx.impactx_pybind.elements.Source
    | impactx.impactx_pybind.elements.SpinMap
    | impactx.impactx_pybind.elements.TaperedPL
    | impactx.impactx_pybind.elements.ThinDipole,
) -> typing.Any:
    """
    Insert an element every s into an element list.

    Returns a new lattice. Elements that are not split are the same objects as
    in the lattice given; an element that a split falls inside is replaced by two
    new elements covering its halves, which carry its parameters but not a Python
    subclass or attributes.
    """
