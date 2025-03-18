"""
Mixin classes for accelerator lattice elements in ImpactX
"""

from __future__ import annotations

__all__ = ["Alignment", "Named", "PipeAperture", "Thick", "Thin"]

class Alignment:
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

class Named:
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
