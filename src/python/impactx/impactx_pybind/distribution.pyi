"""
Particle beam distributions in ImpactX
"""

from __future__ import annotations

__all__ = [
    "Empty",
    "Gaussian",
    "KVdist",
    "Kurth4D",
    "Kurth6D",
    "Semigaussian",
    "Thermal",
    "Triangle",
    "Waterbag",
]

class Empty:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(self) -> None:
        """
        Sets all values to zero.
        """

class Gaussian:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        lambdaX: float,
        lambdaY: float,
        lambdaT: float,
        lambdaPx: float,
        lambdaPy: float,
        lambdaPt: float,
        muxpx: float = 0.0,
        muypy: float = 0.0,
        mutpt: float = 0.0,
    ) -> None:
        """
        A 6D Gaussian distribution
        """

class KVdist:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        lambdaX: float,
        lambdaY: float,
        lambdaT: float,
        lambdaPx: float,
        lambdaPy: float,
        lambdaPt: float,
        muxpx: float = 0.0,
        muypy: float = 0.0,
        mutpt: float = 0.0,
    ) -> None:
        """
        A K-V distribution transversely + a uniform distribution
        in t + a Gaussian distribution in pt
        """

class Kurth4D:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        lambdaX: float,
        lambdaY: float,
        lambdaT: float,
        lambdaPx: float,
        lambdaPy: float,
        lambdaPt: float,
        muxpx: float = 0.0,
        muypy: float = 0.0,
        mutpt: float = 0.0,
    ) -> None:
        """
        A 4D Kurth distribution transversely + a uniform distribution
        in t + a Gaussian distribution in pt
        """

class Kurth6D:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        lambdaX: float,
        lambdaY: float,
        lambdaT: float,
        lambdaPx: float,
        lambdaPy: float,
        lambdaPt: float,
        muxpx: float = 0.0,
        muypy: float = 0.0,
        mutpt: float = 0.0,
    ) -> None:
        """
        A 6D Kurth distribution

        R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)
        C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)
        """

class Semigaussian:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        lambdaX: float,
        lambdaY: float,
        lambdaT: float,
        lambdaPx: float,
        lambdaPy: float,
        lambdaPt: float,
        muxpx: float = 0.0,
        muypy: float = 0.0,
        mutpt: float = 0.0,
    ) -> None:
        """
        A 6D Semi-Gaussian distribution (uniform in position, Gaussian in momentum).
        """

class Thermal:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        k: float,
        kT: float,
        kT_halo: float,
        normalize: float,
        normalize_halo: float,
        halo: float = 0.0,
    ) -> None:
        """
        A stationary thermal or bithermal distribution

        R. D. Ryne, J. Qiang, and A. Adelmann, in Proc. EPAC2004, pp. 1942-1944 (2004)
        """

class Triangle:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        lambdaX: float,
        lambdaY: float,
        lambdaT: float,
        lambdaPx: float,
        lambdaPy: float,
        lambdaPt: float,
        muxpx: float = 0.0,
        muypy: float = 0.0,
        mutpt: float = 0.0,
    ) -> None:
        """
        A triangle distribution for laser-plasma acceleration related applications.

        A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).
        The transverse distribution is a 4D waterbag.
        """

class Waterbag:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs): ...
    def __init__(
        self,
        lambdaX: float,
        lambdaY: float,
        lambdaT: float,
        lambdaPx: float,
        lambdaPy: float,
        lambdaPt: float,
        muxpx: float = 0.0,
        muypy: float = 0.0,
        mutpt: float = 0.0,
    ) -> None:
        """
        A 6D Waterbag distribution
        """
