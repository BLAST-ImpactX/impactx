from __future__ import annotations

import math as math
import typing
import warnings as warnings

import impactx.impactx_pybind
import impactx.impactx_pybind.elements
from amrex.space3d.amrex_3d_pybind import SmallMatrix_6x6_F_SI1_double as Map6x6
from impactx.element_models import select_model, tier_rank, validate_model
from impactx.impactx_pybind import RefPart, elements
from impactx.MADXParser import MADXInputError, MADXParser

__all__: list[str] = [
    "APERTURE_SHAPE_MAP",
    "DRIFT_MODEL_CLASSES",
    "MADXImpactXTranslatorWarning",
    "MADXInputError",
    "MADXParser",
    "Map6x6",
    "QUAD_MODEL_CLASSES",
    "RefPart",
    "beam",
    "elements",
    "lattice",
    "math",
    "read_beam",
    "read_lattice",
    "sc",
    "select_model",
    "tier_rank",
    "validate_model",
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

def _beam_summary(madx_beam) -> str:
    """
    Describe a MAD-X beam for an error or warning.
    """

def _require_finite_attributes(element: dict):
    """
    Stop on a NaN or infinite attribute, as of `L=SQRT(-1)`, of an element in use.

    MAD-X carries it along, but the translation below would treat it like zero.
    """

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
    parsed_beamline,
    nslice=1,
    freq0=0.0,
    options=None,
    ref_mass_MeV=None,
    bv=1.0,
    min_model="linear",
    ref_beta_gamma=None,
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
    :param min_model: lowest ImpactX element model tier to translate into
        ("linear", "paraxial" or "exact"). This raises the lower gate of the
        model choice only. The translator still picks the cheapest model that
        represents the MAD-X element, and a feature that already demands a
        higher tier, such as a thick octupole, keeps it. See `build` below.
    :param ref_beta_gamma: beta*gamma of the reference particle, from the MAD-X
        BEAM total energy. Only needed for models whose strength is expressed
        per unit rigidity rather than in the MAD-X k convention, currently the
        paraxial solenoid. Those models are skipped when this is unknown.
    :return: list of translated dictionaries

    Maintainer note (one-place-of-truth philosophy):
    - Keep element-type mapping in one place (the `d["type"]` chain below).
    - Keep the model-tier choice in one place per element family. Route every
      construction of a tiered element through `build(...)` with a table of the
      tiers ImpactX implements for it, so `min_model` applies uniformly.
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

def read_lattice(madx_file, nslice=1, *, line=None, sequence=None, min_model="linear"):
    """
    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param madx_file: file name to MAD-X file with beamline elements
    :param nslice: number of ds slices per element
    :param line: explicit MAD-X LINE name to expand when no USE command is present
    :param sequence: explicit MAD-X SEQUENCE/PERIOD name to expand
    :param min_model: lowest ImpactX element model tier to translate into
        ("linear", "paraxial" or "exact"), see :py:func:`lattice`
    :return: list of ImpactX.KnownElements
    """

APERTURE_SHAPE_MAP: dict = {
    "rectangle": "rectangular",
    "rectangular": "rectangular",
    "ellipse": "elliptical",
    "circle": "elliptical",
}
DRIFT_MODEL_CLASSES: dict = {
    "linear": impactx.impactx_pybind.elements.Drift,
    "paraxial": impactx.impactx_pybind.elements.ChrDrift,
    "exact": impactx.impactx_pybind.elements.ExactDrift,
}
QUAD_MODEL_CLASSES: dict
