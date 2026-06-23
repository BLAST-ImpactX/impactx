from __future__ import annotations

import math as math
import typing
import warnings as warnings

import impactx.impactx_pybind
from amrex.space3d.amrex_3d_pybind import SmallMatrix_6x6_F_SI1_double as Map6x6
from impactx.impactx_pybind import RefPart, elements
from impactx.MADXParser import MADXParser

__all__: list[str] = [
    "APERTURE_SHAPE_MAP",
    "MADXImpactXTranslatorWarning",
    "MADXParser",
    "Map6x6",
    "RefPart",
    "beam",
    "elements",
    "lattice",
    "math",
    "read_beam",
    "read_lattice",
    "sc",
    "warnings",
]

class MADXImpactXTranslatorWarning(UserWarning):
    """
    Warning category for MAD-X to ImpactX translation fallbacks.
    """

class sc:
    """
    This class is used in lieu of scipy.constants
    to avoid a direct dependency on it.
    At the time of writing, this file was the only
    one requiring scipy in the ImpactX source outside
    of examples.
    """

    c: typing.ClassVar[float] = 299792458.0
    electron_volt: typing.ClassVar[float] = 1.602176634e-19
    m_e: typing.ClassVar[float] = 9.1093837015e-31
    m_p: typing.ClassVar[float] = 1.67262192369e-27
    m_u: typing.ClassVar[float] = 1.6605390666e-27
    physical_constants: typing.ClassVar[dict] = {
        "electron-muon mass ratio": (0.00483633169, "", 1.1e-10)
    }

def _warn(message: str):
    """
    Emit a translation warning without mutating process-wide warning formatting.
    """

def beam(particle, charge=None, mass=None, energy=None):
    """
    Function that converts a list of beam parameter dictionaries in the MADXParser format into ImpactX format

    Rules following https://mad.web.cern.ch/mad/releases/5.02.08/madxuguide.pdf pages 55f.

    :param str particle: reference particle name
    :param float charge: particle charge (proton charge units)
    :param float mass: particle mass in GeV
    :param float energy: total particle energy (GeV)
        - MAD-X default: 1 GeV
    :return dict: dictionary containing particle and beam attributes in ImpactX units
    """

def lattice(
    parsed_beamline, nslice=1, freq0=0.0, options=None, ref_mass_MeV=None, bv=1.0
):
    """
    Function that converts a list of elements in the MADXParser format into ImpactX format
    :param parsed_beamline: list of dictionaries
    :param nslice: number of ds slices per element
    :param freq0: revolution frequency in Hz (from BEAM command)
    :param ref_mass_MeV: reference particle rest mass in MeV (from BEAM command)
    :param bv: MAD-X BEAM `bv` flag (+1 default, -1 to reverse magnet polarities).
        MAD-X exposes this as the per-node ``other_bv`` and multiplies it into
        bend angles, body multipole strengths, solenoid KS/KSI, kicker kicks,
        and RF cavity VOLT (twiss.f90:3856, 4547, 6852; trrun.f90:1610).
    :return: list of translated dictionaries

    Maintainer note (one-place-of-truth philosophy):
    - Keep element-type mapping in one place (`madx_to_impactx_dict` below).
    - Keep bend-body model choice in one place (`make_bend_body_element` below).
    - Prefer composition from existing ImpactX elements over silent drops.
    - Every physics fallback should `_warn(...)` and include a TODO comment.
    """

def read_beam(ref: impactx.impactx_pybind.RefPart, madx_file):
    """
    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements

    .. warning::

       Our MAD-X parser is under active development and provided
       as a preview. Please check any loaded MAD-X beams very
       carefully. Please report your experience and bugs on our
       `issue tracker <https://github.com/BLAST-ImpactX/impactx/issues>`__.

    :param RefPart ref: ImpactX reference particle (passed by reference)
    :param madx_file: file name to MAD-X file with beamline elements
    :return: list of ImpactX.KnownElements
    """

def read_lattice(madx_file, nslice=1, *, line=None, sequence=None):
    """
    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param madx_file: file name to MAD-X file with beamline elements
    :param nslice: number of ds slices per element
    :param line: explicit MAD-X LINE name to expand when no USE command is present
    :param sequence: explicit MAD-X SEQUENCE/PERIOD name to expand
    :return: list of ImpactX.KnownElements
    """

APERTURE_SHAPE_MAP: dict = {
    "rectangle": "rectangular",
    "rectangular": "rectangular",
    "ellipse": "elliptical",
    "circle": "elliptical",
}
