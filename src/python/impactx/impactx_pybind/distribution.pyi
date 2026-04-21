"""
Particle beam distributions in ImpactX
"""

from __future__ import annotations

import typing

__all__: list[str] = [
    "Empty",
    "Gaussian",
    "KVdist",
    "Kurth4D",
    "Kurth6D",
    "Semigaussian",
    "SpinvMF",
    "Thermal",
    "Triangle",
    "Waterbag",
]

class Gaussian:
    def __init__(
        self,
        lambdaX: typing.SupportsFloat | typing.SupportsIndex,
        lambdaY: typing.SupportsFloat | typing.SupportsIndex,
        lambdaT: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPx: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPy: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPt: typing.SupportsFloat | typing.SupportsIndex,
        muxpx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        muypy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        mutpt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        cutX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        cutY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        cutT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A 6D Gaussian distribution
        """

class Kurth4D:
    def __init__(
        self,
        lambdaX: typing.SupportsFloat | typing.SupportsIndex,
        lambdaY: typing.SupportsFloat | typing.SupportsIndex,
        lambdaT: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPx: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPy: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPt: typing.SupportsFloat | typing.SupportsIndex,
        muxpx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        muypy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        mutpt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A 4D Kurth distribution transversely + a uniform distribution
        in t + a Gaussian distribution in pt
        """

class Kurth6D:
    def __init__(
        self,
        lambdaX: typing.SupportsFloat | typing.SupportsIndex,
        lambdaY: typing.SupportsFloat | typing.SupportsIndex,
        lambdaT: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPx: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPy: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPt: typing.SupportsFloat | typing.SupportsIndex,
        muxpx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        muypy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        mutpt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A 6D Kurth distribution

        R. Kurth, Quarterly of Applied Mathematics vol. 32, pp. 325-329 (1978)
        C. Mitchell, K. Hwang and R. D. Ryne, IPAC2021, WEPAB248 (2021)
        """

class KVdist:
    def __init__(
        self,
        lambdaX: typing.SupportsFloat | typing.SupportsIndex,
        lambdaY: typing.SupportsFloat | typing.SupportsIndex,
        lambdaT: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPx: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPy: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPt: typing.SupportsFloat | typing.SupportsIndex,
        muxpx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        muypy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        mutpt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A K-V distribution transversely + a uniform distribution
        in t + a Gaussian distribution in pt
        """

class Empty:
    def __init__(self) -> None:
        """
        Sets all values to zero.
        """

class Semigaussian:
    def __init__(
        self,
        lambdaX: typing.SupportsFloat | typing.SupportsIndex,
        lambdaY: typing.SupportsFloat | typing.SupportsIndex,
        lambdaT: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPx: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPy: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPt: typing.SupportsFloat | typing.SupportsIndex,
        muxpx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        muypy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        mutpt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A 6D Semi-Gaussian distribution (uniform in position, Gaussian in momentum).
        """

class Thermal:
    def __init__(
        self,
        k: typing.SupportsFloat | typing.SupportsIndex,
        kT: typing.SupportsFloat | typing.SupportsIndex,
        kT_halo: typing.SupportsFloat | typing.SupportsIndex,
        normalize: typing.SupportsFloat | typing.SupportsIndex,
        normalize_halo: typing.SupportsFloat | typing.SupportsIndex,
        halo: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A stationary thermal or bithermal distribution

        R. D. Ryne, J. Qiang, and A. Adelmann, in Proc. EPAC2004, pp. 1942-1944 (2004)
        """

class Triangle:
    def __init__(
        self,
        lambdaX: typing.SupportsFloat | typing.SupportsIndex,
        lambdaY: typing.SupportsFloat | typing.SupportsIndex,
        lambdaT: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPx: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPy: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPt: typing.SupportsFloat | typing.SupportsIndex,
        muxpx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        muypy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        mutpt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A triangle distribution for laser-plasma acceleration related applications.

        A ramped, triangular current profile with a Gaussian energy spread (possibly correlated).
        The transverse distribution is a 4D waterbag.
        """

class Waterbag:
    def __init__(
        self,
        lambdaX: typing.SupportsFloat | typing.SupportsIndex,
        lambdaY: typing.SupportsFloat | typing.SupportsIndex,
        lambdaT: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPx: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPy: typing.SupportsFloat | typing.SupportsIndex,
        lambdaPt: typing.SupportsFloat | typing.SupportsIndex,
        muxpx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        muypy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        mutpt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanT: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        meanPt: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispX: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPx: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispY: typing.SupportsFloat | typing.SupportsIndex = 0.0,
        dispPy: typing.SupportsFloat | typing.SupportsIndex = 0.0,
    ) -> None:
        """
        A 6D Waterbag distribution
        """

class SpinvMF:
    @staticmethod
    def inverse_Langevin(pmag: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
        This function evaluates the inverse Langevin function, in order to return the value of concentration (kappa) required to produce a given polarization magnitude.
        """
    def __init__(
        self,
        mux: typing.SupportsFloat | typing.SupportsIndex,
        muy: typing.SupportsFloat | typing.SupportsIndex,
        muz: typing.SupportsFloat | typing.SupportsIndex,
    ) -> None:
        """
        A von Mises-Fisher (vMF) distribution on the unit 2-sphere, for particle spin.
        """
