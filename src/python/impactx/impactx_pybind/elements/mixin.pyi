"""
Mixin classes for accelerator lattice elements in ImpactX
"""

from __future__ import annotations

import amrex.space3d.amrex_3d_pybind

__all__ = ["Alignment", "LinearTransport", "Named", "PipeAperture", "Thick", "Thin"]

class Alignment:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    @property
    def dx(self) -> float:
        """
        horizontal translation error in m
        """
    @dx.setter
    def dx(self, arg1: float) -> None: ...
    @property
    def dy(self) -> float:
        """
        vertical translation error in m
        """
    @dy.setter
    def dy(self, arg1: float) -> None: ...
    @property
    def rotation(self) -> float:
        """
        rotation error in the transverse plane in degree
        """
    @rotation.setter
    def rotation(self, arg1: float) -> None: ...

class LinearTransport:
    Map6x6 = amrex.space3d.amrex_3d_pybind.SmallMatrix_6x6_F_SI1_double
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...

class Named:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    @property
    def has_name(self) -> bool: ...
    @property
    def name(self) -> str:
        """
        segment length in m
        """
    @name.setter
    def name(self, arg1: str) -> None: ...

class PipeAperture:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    @property
    def aperture_x(self) -> float:
        """
        horizontal aperture in m
        """
    @property
    def aperture_y(self) -> float:
        """
        vertical aperture in m
        """

class Thick:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    @property
    def ds(self) -> float:
        """
        segment length in m
        """
    @ds.setter
    def ds(self, arg1: float) -> None: ...
    @property
    def nslice(self) -> int:
        """
        number of slices used for the application of space charge
        """
    @nslice.setter
    def nslice(self, arg1: int) -> None: ...

class Thin:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    @property
    def ds(self) -> float:
        """
        segment length in m
        """
    @property
    def nslice(self) -> int:
        """
        number of slices used for the application of space charge
        """
