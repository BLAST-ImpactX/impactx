"""
Accelerator lattice elements in ImpactX
"""

from __future__ import annotations

import collections.abc
import typing

import amrex.space3d.amrex_3d_pybind
import impactx.impactx_pybind

from . import mixin, transformation

__all__: list[str] = [
    "Aperture",
    "BeamMonitor",
    "Buncher",
    "CFbend",
    "ChrAcc",
    "ChrDrift",
    "ChrPlasmaLens",
    "ChrQuad",
    "ConstF",
    "DipEdge",
    "Drift",
    "Empty",
    "ExactCFbend",
    "ExactDrift",
    "ExactMultipole",
    "ExactQuad",
    "ExactSbend",
    "Kicker",
    "KnownElementsList",
    "LinearMap",
    "Marker",
    "Multipole",
    "NonlinearLens",
    "PRot",
    "PlaneXYRot",
    "Programmable",
    "Quad",
    "QuadEdge",
    "RFCavity",
    "Sbend",
    "ShortRF",
    "SoftQuadrupole",
    "SoftSolenoid",
    "Sol",
    "Source",
    "TaperedPL",
    "ThinDipole",
    "mixin",
    "transformation",
]

class Aperture(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        aperture_x: typing.SupportsFloat,
        aperture_y: typing.SupportsFloat,
        repeat_x: typing.SupportsFloat = 0,
        repeat_y: typing.SupportsFloat = 0,
        shift_odd_x: bool = False,
        shape: str = "rectangular",
        action: str = "transmit",
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A short collimator element applying a transverse aperture boundary.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def action(self) -> str:
        """
        action type (transmit, absorb)
        """
    @action.setter
    def action(self, arg1: str) -> None: ...
    @property
    def aperture_x(self) -> float:
        """
        maximum horizontal coordinate
        """
    @aperture_x.setter
    def aperture_x(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def aperture_y(self) -> float:
        """
        maximum vertical coordinate
        """
    @aperture_y.setter
    def aperture_y(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def repeat_x(self) -> float:
        """
        horizontal period for repeated aperture masking
        """
    @repeat_x.setter
    def repeat_x(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def repeat_y(self) -> float:
        """
        vertical period for repeated aperture masking
        """
    @repeat_y.setter
    def repeat_y(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def shape(self) -> str:
        """
        aperture type (rectangular, elliptical)
        """
    @shape.setter
    def shape(self, arg1: str) -> None: ...
    @property
    def shift_odd_x(self) -> bool:
        """
        for hexagonal/triangular mask patterns: horizontal shift of every 2nd (odd) vertical period by repeat_x / 2. Use alignment offsets dx,dy to move whole mask as needed.
        """
    @shift_odd_x.setter
    def shift_odd_x(self, arg1: bool) -> None: ...

class BeamMonitor(mixin.Thin):
    def __init__(
        self,
        name: str,
        backend: str = "default",
        encoding: str = "g",
        period_sample_intervals: typing.SupportsInt = 1,
    ) -> None:
        """
        This element writes the particle beam out to openPMD data.
        """
    def __repr__(self) -> str: ...
    def finalize(self) -> None: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def alpha(self) -> float:
        """
        Twiss alpha of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.
        """
    @alpha.setter
    def alpha(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def beta(self) -> float:
        """
        Twiss beta (in meters) of the bare linear lattice at the location of output for the nonlinear IOTA invariants H and I.
        Horizontal and vertical values must be equal.
        """
    @beta.setter
    def beta(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def cn(self) -> float:
        """
        Scale factor (in meters^(1/2)) of the IOTA nonlinear magnetic insert element used for computing H and I.
        """
    @cn.setter
    def cn(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def name(self) -> str:
        """
        name of the series
        """
    @property
    def nonlinear_lens_invariants(self) -> bool:
        """
        Compute and output the invariants H and I within the nonlinear magnetic insert element
        """
    @nonlinear_lens_invariants.setter
    def nonlinear_lens_invariants(self, arg1: bool) -> None: ...
    @property
    def tn(self) -> float:
        """
        Dimensionless strength of the IOTA nonlinear magnetic insert element used for computing H and I.
        """
    @tn.setter
    def tn(self, arg1: typing.SupportsFloat) -> None: ...

class Buncher(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        V: typing.SupportsFloat,
        k: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A short linear RF cavity element at zero-crossing for bunching.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def V(self) -> float:
        """
        Normalized RF voltage drop V = Emax*L/(c*Brho)
        """
    @V.setter
    def V(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def k(self) -> float:
        """
        Wavenumber of RF in 1/m
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...

class CFbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        rc: typing.SupportsFloat,
        k: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal combined function bend (sector bend with quadrupole component).
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def k(self) -> float:
        """
        Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def rc(self) -> float:
        """
        Radius of curvature in m
        """
    @rc.setter
    def rc(self, arg1: typing.SupportsFloat) -> None: ...

class ChrAcc(mixin.Named, mixin.Thick, mixin.Alignment):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        ez: typing.SupportsFloat,
        bz: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A region of Uniform Acceleration, with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def bz(self) -> float:
        """
        magnetic field strength in 1/m
        """
    @bz.setter
    def bz(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def ez(self) -> float:
        """
        electric field strength in 1/m
        """
    @ez.setter
    def ez(self, arg1: typing.SupportsFloat) -> None: ...

class ChrDrift(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A Drift with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...

class ChrPlasmaLens(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        k: typing.SupportsFloat,
        unit: typing.SupportsInt = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        An active Plasma Lens with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def k(self) -> float:
        """
        focusing strength in 1/m^2 (or T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for focusing strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt) -> None: ...

class ChrQuad(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        k: typing.SupportsFloat,
        unit: typing.SupportsInt = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A Quadrupole magnet with chromatic effects included.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def k(self) -> float:
        """
        quadrupole strength in 1/m^2 (or T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for quad strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt) -> None: ...

class ConstF(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        kx: typing.SupportsFloat,
        ky: typing.SupportsFloat,
        kt: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A linear Constant Focusing element.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def kt(self) -> float:
        """
        focusing t strength in 1/m
        """
    @kt.setter
    def kt(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def kx(self) -> float:
        """
        focusing x strength in 1/m
        """
    @kx.setter
    def kx(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def ky(self) -> float:
        """
        focusing y strength in 1/m
        """
    @ky.setter
    def ky(self, arg1: typing.SupportsFloat) -> None: ...

class DipEdge(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        psi: typing.SupportsFloat,
        rc: typing.SupportsFloat,
        g: typing.SupportsFloat,
        K2: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        Edge focusing associated with bend entry or exit.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def K2(self) -> float:
        """
        Fringe field integral (unitless)
        """
    @K2.setter
    def K2(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def g(self) -> float:
        """
        Gap parameter in m
        """
    @g.setter
    def g(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def psi(self) -> float:
        """
        Pole face angle in rad
        """
    @psi.setter
    def psi(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def rc(self) -> float:
        """
        Radius of curvature in m
        """
    @rc.setter
    def rc(self, arg1: typing.SupportsFloat) -> None: ...

class Drift(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A drift.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...

class Empty(mixin.Thin):
    def __init__(self) -> None:
        """
        This element does nothing.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...

class ExactCFbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        k_normal: collections.abc.Sequence[typing.SupportsFloat],
        k_skew: collections.abc.Sequence[typing.SupportsFloat],
        unit: typing.SupportsInt = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        int_order: typing.SupportsInt = 2,
        mapsteps: typing.SupportsInt = 5,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A thick combined function bending magnet using the exact nonlinear Hamiltonian.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def int_order(self) -> int:
        """
        order of symplectic integration used for particle push in applied fields
        """
    @int_order.setter
    def int_order(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for particle push in the applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for multipole strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt) -> None: ...

class ExactDrift(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A Drift using the exact nonlinear map.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...

class ExactMultipole(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        k_normal: collections.abc.Sequence[typing.SupportsFloat],
        k_skew: collections.abc.Sequence[typing.SupportsFloat],
        unit: typing.SupportsInt = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        int_order: typing.SupportsInt = 2,
        mapsteps: typing.SupportsInt = 5,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A thick Multipole magnet using the exact nonlinear Hamiltonian.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def int_order(self) -> int:
        """
        order of symplectic integration used for particle push in applied fields
        """
    @int_order.setter
    def int_order(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for particle push in the applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for multipole strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt) -> None: ...

class ExactQuad(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        k: typing.SupportsFloat,
        unit: typing.SupportsInt = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        int_order: typing.SupportsInt = 2,
        mapsteps: typing.SupportsInt = 5,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A Quadrupole magnet using the exact nonlinear Hamiltonian.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def int_order(self) -> int:
        """
        order of symplectic integration used for particle push in applied fields
        """
    @int_order.setter
    def int_order(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def k(self) -> float:
        """
        quadrupole strength in 1/m^2 (or T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for particle push in the applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for quad strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt) -> None: ...

class ExactSbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        phi: typing.SupportsFloat,
        B: typing.SupportsFloat = 0.0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal sector bend using the exact nonlinear map.  When B = 0, the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def rc(self, ref: impactx.impactx_pybind.RefPart) -> float:
        """
        Radius of curvature in m
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def B(self) -> float:
        """
        Magnetic field in Tesla; when B = 0 (default), the reference bending radius is defined by r0 = length / (angle in rad), corresponding to a magnetic field of B = rigidity / r0; otherwise the reference bending radius is defined by r0 = rigidity / B
        """
    @B.setter
    def B(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def phi(self) -> float:
        """
        Bend angle in degrees
        """
    @phi.setter
    def phi(self, arg1: typing.SupportsFloat) -> None: ...

class Kicker(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        xkick: typing.SupportsFloat,
        ykick: typing.SupportsFloat,
        unit: str = "dimensionless",
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A thin transverse kicker element. Kicks are for unit "dimensionless" or in "T-m".
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def xkick(self) -> float:
        """
        horizontal kick strength (dimensionless OR T-m)
        """
    @xkick.setter
    def xkick(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def ykick(self) -> float:
        """
        vertical kick strength (dimensionless OR T-m)
        """
    @ykick.setter
    def ykick(self, arg1: typing.SupportsFloat) -> None: ...

class KnownElementsList:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: impactx.impactx_pybind.elements.Empty
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
        | impactx.impactx_pybind.elements.TaperedPL
        | impactx.impactx_pybind.elements.ThinDipole,
    ) -> None: ...
    @typing.overload
    def __init__(self, arg0: list) -> None: ...
    def __iter__(
        self,
    ) -> collections.abc.Iterator[
        impactx.impactx_pybind.elements.Empty
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
        | impactx.impactx_pybind.elements.TaperedPL
        | impactx.impactx_pybind.elements.ThinDipole
    ]: ...
    def __len__(self) -> int:
        """
        The length of the list.
        """
    def append(
        self,
        arg0: impactx.impactx_pybind.elements.Empty
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
        | impactx.impactx_pybind.elements.TaperedPL
        | impactx.impactx_pybind.elements.ThinDipole,
    ) -> None:
        """
        Add a single element to the list.
        """
    def clear(self) -> None:
        """
        Clear the list to become empty.
        """
    @typing.overload
    def extend(self, arg0: KnownElementsList) -> KnownElementsList:
        """
        Add a list of elements to the list.
        """
    @typing.overload
    def extend(self, arg0: list) -> KnownElementsList:
        """
        Add a list of elements to the list.
        """
    def from_pals(self, pals_beamline, nslice=1):
        """
        Load and append a lattice from a Particle Accelerator Lattice Standard (PALS) Python BeamLine.

            https://github.com/campa-consortium/pals-python

        """
    def load_file(self, filename, nslice=1):
        """
        Load and append a lattice file from MAD-X (.madx) or PALS (e.g., .pals.yaml) formats.
        """
    def plot_survey(
        self, ref=None, ax=None, legend=True, legend_ncols=5, palette="cern-lhc"
    ):
        """
        Plot over s of all elements in the KnownElementsList.

            A positive element strength denotes horizontal focusing (e.g. for quadrupoles) and bending to the right (for dipoles).  In general, this depends on both the sign of the field and the sign of the charge.

            Parameters
            ----------
            self : ImpactXParticleContainer_*
                The KnownElementsList class in ImpactX
            ref : RefPart
                A reference particle, checked for the charge sign to plot focusing/defocusing strength directions properly.
            ax : matplotlib axes
                A plotting area in matplotlib (called axes there).
            legend: bool
                Plot a legend if true.
            legend_ncols: int
                Number of columns for lattice element types in the legend.
            palette: string
                Color palette.

            Returns
            -------
            Either populates the matplotlib axes in ax or creates a new axes containing the plot.

        """
    def pop_back(self) -> None:
        """
        Return and remove the last element of the list.
        """

class LinearMap(mixin.Named, mixin.Alignment):
    def __init__(
        self,
        R: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ds: typing.SupportsFloat = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        (A user-provided linear map, represented as a 6x6 transport matrix.)
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def R(self) -> amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double:
        """
        linear map as a 6x6 transport matrix
        """
    @R.setter
    def R(
        self, arg1: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
    ) -> None: ...
    @property
    def ds(self) -> float:
        """
        segment length in m
        """
    @ds.setter
    def ds(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def nslice(self) -> int:
        """
        one, because we do not support slicing of this element
        """

class Marker(mixin.Named, mixin.Thin):
    def __init__(self, name: str) -> None:
        """
        This named element does nothing.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...

class Multipole(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        multipole: typing.SupportsInt,
        K_normal: typing.SupportsFloat,
        K_skew: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A general thin multipole element.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def K_normal(self) -> float:
        """
        Integrated normal multipole coefficient (1/meter^m)
        """
    @K_normal.setter
    def K_normal(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def K_skew(self) -> float:
        """
        Integrated skew multipole coefficient (1/meter^m)
        """
    @K_skew.setter
    def K_skew(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def multipole(self) -> int:
        """
        index m (m=1 dipole, m=2 quadrupole, m=3 sextupole etc.)
        """
    @multipole.setter
    def multipole(self, arg1: typing.SupportsFloat) -> None: ...

class NonlinearLens(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        knll: typing.SupportsFloat,
        cnll: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        Single short segment of the nonlinear magnetic insert element.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def cnll(self) -> float:
        """
        distance of singularities from the origin (m)
        """
    @cnll.setter
    def cnll(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def knll(self) -> float:
        """
        integrated strength of the nonlinear lens (m)
        """
    @knll.setter
    def knll(self, arg1: typing.SupportsFloat) -> None: ...

class PRot(mixin.Named, mixin.Thin):
    def __init__(
        self,
        phi_in: typing.SupportsFloat,
        phi_out: typing.SupportsFloat,
        name: str | None = None,
    ) -> None:
        """
        An exact pole-face rotation in the x-z plane. Both angles are in degrees.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def phi_in(self) -> float:
        """
        angle of the reference particle with respect to the longitudinal (z) axis in the original frame in degrees
        """
    @phi_in.setter
    def phi_in(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def phi_out(self) -> float:
        """
        angle of the reference particle with respect to the longitudinal (z) axis in the rotated frame in degrees
        """
    @phi_out.setter
    def phi_out(self, arg1: typing.SupportsFloat) -> None: ...

class PlaneXYRot(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        angle: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A rotation in the x-y plane.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def angle(self) -> float:
        """
        Rotation angle (rad).
        """
    @angle.setter
    def angle(self, arg1: typing.SupportsFloat) -> None: ...

class Programmable(mixin.Named):
    def __init__(
        self,
        ds: typing.SupportsFloat = 0.0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A programmable beam optics element.
        """
    def __repr__(self) -> str: ...
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def beam_particles(
        self,
    ) -> collections.abc.Callable[
        [impactx.impactx_pybind.ImpactXParIter, impactx.impactx_pybind.RefPart], None
    ]:
        """
        hook for beam particles (pti, RefPart)
        """
    @beam_particles.setter
    def beam_particles(
        self,
        arg1: collections.abc.Callable[
            [impactx.impactx_pybind.ImpactXParIter, impactx.impactx_pybind.RefPart],
            None,
        ],
    ) -> None: ...
    @property
    def ds(self) -> float: ...
    @ds.setter
    def ds(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def nslice(self) -> int: ...
    @nslice.setter
    def nslice(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def push(
        self,
    ) -> collections.abc.Callable[
        [
            impactx.impactx_pybind.ImpactXParticleContainer,
            typing.SupportsInt,
            typing.SupportsInt,
        ],
        None,
    ]:
        """
        hook for push of whole container (pc, step, period)
        """
    @push.setter
    def push(
        self,
        arg1: collections.abc.Callable[
            [
                impactx.impactx_pybind.ImpactXParticleContainer,
                typing.SupportsInt,
                typing.SupportsInt,
            ],
            None,
        ],
    ) -> None: ...
    @property
    def ref_particle(
        self,
    ) -> collections.abc.Callable[[impactx.impactx_pybind.RefPart], None]:
        """
        hook for reference particle (RefPart)
        """
    @ref_particle.setter
    def ref_particle(
        self, arg1: collections.abc.Callable[[impactx.impactx_pybind.RefPart], None]
    ) -> None: ...
    @property
    def threadsafe(self) -> bool:
        """
        allow threading via OpenMP for the particle iterator loop, default=False (note: if OMP backend is active)
        """
    @threadsafe.setter
    def threadsafe(self, arg1: bool) -> None: ...

class Quad(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        k: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A Quadrupole magnet.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def k(self) -> float:
        """
        Quadrupole strength in m^(-2) (MADX convention) = (gradient in T/m) / (rigidity in T-m) k > 0 horizontal focusing k < 0 horizontal defocusing
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...

class QuadEdge(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        k: typing.SupportsFloat,
        unit: typing.SupportsInt = 0,
        flag: str = "entry",
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A thin quadrupole fringe field element. Flag must be "entry" or "exit".
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def k(self) -> float:
        """
        quadrupole focusing strength (1/meter^2 OR T/m)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def unit(self) -> int:
        """
        unit specification for quad strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt) -> None: ...

class RFCavity(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        escale: typing.SupportsFloat,
        freq: typing.SupportsFloat,
        phase: typing.SupportsFloat,
        cos_coefficients: collections.abc.Sequence[typing.SupportsFloat],
        sin_coefficients: collections.abc.Sequence[typing.SupportsFloat],
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        mapsteps: typing.SupportsInt = 1,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        An RF cavity.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def escale(self) -> float:
        """
        scaling factor for on-axis RF electric field in 1/m = (peak on-axis electric field Ez in MV/m) / (particle rest energy in MeV)
        """
    @escale.setter
    def escale(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def freq(self) -> float:
        """
        RF frequency in Hz
        """
    @freq.setter
    def freq(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def phase(self) -> float:
        """
        RF driven phase in degrees
        """
    @phase.setter
    def phase(self, arg1: typing.SupportsFloat) -> None: ...

class Sbend(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        rc: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal sector bend.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def rc(self, ref: impactx.impactx_pybind.RefPart = None) -> float:
        """
        Radius of curvature in m
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...

class ShortRF(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        V: typing.SupportsFloat,
        freq: typing.SupportsFloat,
        phase: typing.SupportsFloat = -90.0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A short RF cavity element.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def V(self) -> float:
        """
        Normalized RF voltage V = maximum energy gain/(m*c^2)
        """
    @V.setter
    def V(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def freq(self) -> float:
        """
        RF frequency in Hz
        """
    @freq.setter
    def freq(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def phase(self) -> float:
        """
        RF synchronous phase in degrees (phase = 0 corresponds to maximum energy gain, phase = -90 corresponds go zero energy gain for bunching)
        """
    @phase.setter
    def phase(self, arg1: typing.SupportsFloat) -> None: ...

class SoftQuadrupole(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        gscale: typing.SupportsFloat,
        cos_coefficients: collections.abc.Sequence[typing.SupportsFloat],
        sin_coefficients: collections.abc.Sequence[typing.SupportsFloat],
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        mapsteps: typing.SupportsInt = 1,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A soft-edge quadrupole.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def gscale(self) -> float:
        """
        Scaling factor for on-axis field gradient in inverse meters
        """
    @gscale.setter
    def gscale(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt) -> None: ...

class SoftSolenoid(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        bscale: typing.SupportsFloat,
        cos_coefficients: collections.abc.Sequence[typing.SupportsFloat],
        sin_coefficients: collections.abc.Sequence[typing.SupportsFloat],
        unit: typing.SupportsFloat = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        mapsteps: typing.SupportsInt = 1,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        A soft-edge solenoid.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def bscale(self) -> float:
        """
        Scaling factor for on-axis magnetic field Bz in inverse meters (if unit = 0) or magnetic field Bz in T (SI units, if unit = 1)
        """
    @bscale.setter
    def bscale(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def mapsteps(self) -> int:
        """
        number of integration steps per slice used for map and reference particle push in applied fields
        """
    @mapsteps.setter
    def mapsteps(self, arg1: typing.SupportsInt) -> None: ...
    @property
    def unit(self) -> int:
        """
        specification of units for scaling of the on-axis longitudinal magnetic field
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsFloat) -> None: ...

class Sol(mixin.Named, mixin.Thick, mixin.Alignment, mixin.PipeAperture):
    def __init__(
        self,
        ds: typing.SupportsFloat,
        ks: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        aperture_x: typing.SupportsFloat = 0,
        aperture_y: typing.SupportsFloat = 0,
        nslice: typing.SupportsInt = 1,
        name: str | None = None,
    ) -> None:
        """
        An ideal hard-edge Solenoid magnet.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def ks(self) -> float:
        """
        Solenoid strength in m^(-1) (MADX convention) in (magnetic field Bz in T) / (rigidity in T-m)
        """
    @ks.setter
    def ks(self, arg1: typing.SupportsFloat) -> None: ...

class Source(mixin.Named, mixin.Thin):
    def __init__(
        self, distribution: str, openpmd_path: str, name: str | None = None
    ) -> None:
        """
        A particle source.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def distribution(self) -> str:
        """
        Distribution type of particles in the source
        """
    @distribution.setter
    def distribution(self, arg1: str) -> None: ...
    @property
    def series_name(self) -> str:
        """
        Path to openPMD series as accepted by openPMD_api.Series
        """
    @series_name.setter
    def series_name(self, arg1: str) -> None: ...

class TaperedPL(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        k: typing.SupportsFloat,
        taper: typing.SupportsFloat,
        unit: typing.SupportsInt = 0,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A thin nonlinear plasma lens with transverse (horizontal) taper

                     .. math::

                        B_x = g \\left( y + \\frac{xy}{D_x} \\right), \\quad \\quad B_y = -g \\left(x + \\frac{x^2 + y^2}{2 D_x} \\right)

                     where :math:`g` is the (linear) field gradient in T/m and :math:`D_x` is the targeted horizontal dispersion in m.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def k(self) -> float:
        """
        integrated focusing strength in m^(-1) (if unit = 0) or integrated focusing strength in T (if unit = 1)
        """
    @k.setter
    def k(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def taper(self) -> float:
        """
        horizontal taper parameter in m^(-1) = 1 / (target horizontal dispersion in m)
        """
    @taper.setter
    def taper(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def unit(self) -> int:
        """
        specification of units for plasma lens focusing strength
        """
    @unit.setter
    def unit(self, arg1: typing.SupportsInt) -> None: ...

class ThinDipole(mixin.Named, mixin.Thin, mixin.Alignment):
    def __init__(
        self,
        theta: typing.SupportsFloat,
        rc: typing.SupportsFloat,
        dx: typing.SupportsFloat = 0,
        dy: typing.SupportsFloat = 0,
        rotation: typing.SupportsFloat = 0,
        name: str | None = None,
    ) -> None:
        """
        A thin kick model of a dipole bend.
        """
    def __repr__(self) -> str: ...
    @typing.overload
    def push(
        self,
        pc: impactx.impactx_pybind.ImpactXParticleContainer,
        step: typing.SupportsInt = 0,
        period: typing.SupportsInt = 0,
    ) -> None:
        """
        Push first the reference particle, then all other particles.
        """
    @typing.overload
    def push(
        self,
        cm: amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double,
        ref: impactx.impactx_pybind.RefPart,
    ) -> None:
        """
        Linear push of the covariance matrix through an element. Expects that the reference particle was advanced first.
        """
    def to_dict(
        self,
    ) -> dict[
        str,
        float
        | int
        | int
        | str
        | list[float]
        | list[int]
        | list[int]
        | amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
        | None,
    ]: ...
    @property
    def rc(self) -> float:
        """
        Effective curvature radius (meters)
        """
    @rc.setter
    def rc(self, arg1: typing.SupportsFloat) -> None: ...
    @property
    def theta(self) -> float:
        """
        Bend angle (degrees)
        """
    @theta.setter
    def theta(self, arg1: typing.SupportsFloat) -> None: ...
