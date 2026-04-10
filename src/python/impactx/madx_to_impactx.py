#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math
import warnings

from impactx import RefPart, elements

from .MADXParser import MADXParser

WARN_PREFIX = "[WARNING MADX-ImpactX-Translator] "
TRACKED_TRANSLATION_OPTIONS = ("rbarc", "thin_foc")
# Centralized aperture shape mapping (MAD-X -> ImpactX Aperture.shape).
# Edit this table first when adding support for new MAD-X aperture types.
APERTURE_SHAPE_MAP = {
    "rectangle": "rectangular",
    "rectangular": "rectangular",
    "lhcscreen": "rectangular",
    "ellipse": "elliptical",
    "circle": "elliptical",
}


class MADXImpactXTranslatorWarning(UserWarning):
    """Warning category for MAD-X to ImpactX translation fallbacks."""

    pass


# Custom format for translator warnings: omit Python source-line echo.
_translator_original_formatwarning = warnings.formatwarning


def _translator_formatwarning(message, category, filename, lineno, line=None):
    if issubclass(category, MADXImpactXTranslatorWarning):
        return f"{WARN_PREFIX}{str(message)}\n"
    return _translator_original_formatwarning(message, category, filename, lineno, line)


warnings.formatwarning = _translator_formatwarning


def _warn(message: str):
    """Emit a standardized translation warning."""
    warnings.warn(message, MADXImpactXTranslatorWarning)


class sc:
    """
    This class is used in lieu of scipy.constants
    to avoid a direct dependency on it.
    At the time of writing, this file was the only
    one requiring scipy in the ImpactX source outside
    of examples.
    """

    c = 299792458.0
    electron_volt = 1.602176634e-19
    physical_constants = {"electron-muon mass ratio": (0.00483633169, "", 1.1e-10)}
    m_e = 9.1093837015e-31
    m_p = 1.67262192369e-27
    m_u = 1.6605390666e-27


def lattice(parsed_beamline, nslice=1, freq0=0.0, options=None):
    """
    Function that converts a list of elements in the MADXParser format into ImpactX format
    :param parsed_beamline: list of dictionaries
    :param nslice: number of ds slices per element
    :param freq0: revolution frequency in MHz (from BEAM command)
    :return: list of translated dictionaries

    Maintainer note (one-place-of-truth philosophy):
    - Keep element-type mapping in one place (`madx_to_impactx_dict` below).
    - Keep bend-body model choice in one place (`make_bend_body_element` below).
    - Prefer composition from existing ImpactX elements over silent drops.
    - Every physics fallback should `_warn(...)` and include a TODO comment.
    """

    # mapping is "MAD-X": "ImpactX",
    madx_to_impactx_dict = {
        "MARKER": "None",
        "BEAMBEAM": "None",  # beam-beam interaction, no ImpactX equivalent
        "DRIFT": "Drift",
        "SBEND": "Sbend",  # Sector Bending Magnet
        "RBEND": "Sbend",  # Rectangular Bending Magnet -> DipEdge + Sbend + DipEdge
        "SOLENOID": "Sol",  # Ideal, thick Solenoid: MAD-X user guide 10.9 p78
        "QUADRUPOLE": "Quad",  # Quadrupole
        "DIPEDGE": "DipEdge",
        # Kicker, idealized thin element,
        # MADX also defines length "L" and a roll angle around the longitudinal axis "TILT"
        # https://mad.web.cern.ch/mad/webguide/manual.html#Ch11.S11
        "KICKER": "Kicker",
        "TKICKER": "Kicker",  # thin kicker, treated same as KICKER
        "HKICKER": "Kicker",  # horizontal kicker
        "VKICKER": "Kicker",  # vertical kicker
        # note: in MAD-X, this keeps track only of the beam centroid,
        # "In addition it serves to record the beam position for closed orbit correction."
        "MONITOR": "BeamMonitor",  # drift + output diagnostics
        "HMONITOR": "BeamMonitor",  # horizontal monitor
        "VMONITOR": "BeamMonitor",  # vertical monitor
        "MULTIPOLE": "Multipole",
        "SEXTUPOLE": "ExactMultipole",
        "OCTUPOLE": "ExactMultipole",
        "RFCAVITY": "RFCavity",
        "NLLENS": "NonlinearLens",
        "COLLIMATOR": "Aperture",  # fallback: geometry only (+ drift for length)
        "RCOLLIMATOR": "Aperture",  # rectangular collimator alias
        "ECOLLIMATOR": "Aperture",  # elliptical collimator alias
        "PLACEHOLDER": "Marker",  # non-physical placement/device placeholder
        "INSTRUMENT": "Marker",  # optically drift-like diagnostics placeholder
        # TODO Figure out how to identify these
        "ShortRF": "ShortRF",
        "ConstF": "ConstF",
    }

    impactx_beamline = []
    options = options or {}
    # MAD-X defaults from OPTION command table:
    # rbarc=true, thin_foc=true
    rbarc = bool(options.get("rbarc", True))
    thin_foc = bool(options.get("thin_foc", True))
    rad_to_deg = 180.0 / math.pi
    # MAD-X bend maps are driven by ANGLE/TILT; K0/K0S are obsolete for map construction.
    # H1/H2 (pole-face curvature) still affect MAD-X fringe modeling and are not translated here.
    unsupported_bend_attrs = ("k0", "k0s", "h1", "h2")
    warned_once = set()

    class _TrackedElement(dict):
        """Dictionary wrapper that records accessed keys."""

        def __init__(self, data):
            super().__init__(data)
            self.accessed = set()

        def __getitem__(self, key):
            self.accessed.add(key)
            return super().__getitem__(key)

        def get(self, key, default=None):
            self.accessed.add(key)
            return super().get(key, default)

        def __contains__(self, key):
            self.accessed.add(key)
            return super().__contains__(key)

    def warn_once(key, message):
        """Emit a warning only once per lattice translation call."""
        if key not in warned_once:
            _warn(message)
            warned_once.add(key)

    def warn_unread_element_attributes(elem):
        """Warn if a MAD-X element attribute is present but not consumed."""
        ignored_meta = {"name", "type", "at", "from"}
        etype = elem.get("type", "")
        for attr in elem.keys():
            if attr in ignored_meta:
                continue
            if attr not in elem.accessed:
                warn_once(
                    f"unread_attr:{etype}:{attr}",
                    f"{etype.upper()} attribute '{attr}' is present in MAD-X input but not consumed by translator.",
                )

    def append_thin_with_length(thin_element, length, element_name):
        """Model finite-length thin-lens MAD-X elements as drift-kick-drift."""
        if length > 0.0:
            # TODO(audit): This is a paraxial length-preserving approximation.
            # Replace with dedicated thick elements when available.
            impactx_beamline.append(
                elements.Drift(
                    name=element_name + "_drift_in", ds=0.5 * length, nslice=nslice
                )
            )
            impactx_beamline.append(thin_element)
            impactx_beamline.append(
                elements.Drift(
                    name=element_name + "_drift_out", ds=0.5 * length, nslice=nslice
                )
            )
        else:
            impactx_beamline.append(thin_element)

    def parse_aperture_half_axes(aperture):
        """Extract (aperture_x, aperture_y) half-axes from MAD-X aperture data."""
        ax = 0.0
        ay = 0.0
        if isinstance(aperture, (list, tuple)):
            if len(aperture) >= 1:
                ax = float(aperture[0])
            if len(aperture) >= 2:
                ay = float(aperture[1])
        elif isinstance(aperture, (int, float)):
            ax = float(aperture)
            ay = float(aperture)
        return ax, ay

    def parse_aperture_offsets(aper_offset):
        """Extract (dx, dy) aperture offsets from MAD-X aper_offset data."""
        dx = 0.0
        dy = 0.0
        if isinstance(aper_offset, (list, tuple)):
            if len(aper_offset) >= 1:
                dx = float(aper_offset[0])
            if len(aper_offset) >= 2:
                dy = float(aper_offset[1])
        return dx, dy

    def aperture_params_from_madx(elem):
        """Build aperture parameter dict from MAD-X aperture metadata.

        Returns:
            dict with Aperture constructor kwargs or None if unsupported data.
        """
        apertype = str(elem.get("apertype", "")).lower()
        aperture = elem.get("aperture", [])
        aper_offset = elem.get("aper_offset", [0.0, 0.0])

        if not (apertype or aperture):
            return None

        shape = APERTURE_SHAPE_MAP.get(apertype)
        ax, ay = parse_aperture_half_axes(aperture)
        dx, dy = parse_aperture_offsets(aper_offset)

        # Circle in MAD-X uses one radius; mirror to y if only x is provided.
        if apertype == "circle" and ay == 0.0:
            ay = ax

        if shape is None or ax <= 0.0 or ay <= 0.0:
            return None

        return {
            "name": elem["name"] + "_aperture",
            "aperture_x": ax,
            "aperture_y": ay,
            "shape": shape,
            "dx": dx,
            "dy": dy,
        }

    def aperture_element_from_madx(elem):
        """Build an ImpactX Aperture element from MAD-X aperture metadata."""
        params = aperture_params_from_madx(elem)
        if params is None:
            return None
        return elements.Aperture(**params)

    def marker_aperture_element(marker):
        """Build an ImpactX Aperture element from MARKER aperture metadata."""
        return aperture_element_from_madx(marker)

    def append_collimator_equivalent(*, name, ds, aperture_params):
        """Append a MAD-X collimator equivalent using current ImpactX primitives.

        MAD-X semantics (manual 11.21):
        - Optically behaves like a drift.
        - Aperture checks are applied at entrance (and exit if L>0).

        ImpactX has no dedicated collimator jaw-material model here, so we
        preserve drift optics and approximate boundary checks with Aperture
        elements at both ends for finite-length collimators.
        """
        if ds > 0.0 and aperture_params is not None:
            # Create distinct aperture elements so entry/exit checks are explicit.
            impactx_beamline.append(
                elements.Aperture(
                    name=name + "_aperture_entry",
                    aperture_x=aperture_params["aperture_x"],
                    aperture_y=aperture_params["aperture_y"],
                    shape=aperture_params["shape"],
                    dx=aperture_params["dx"],
                    dy=aperture_params["dy"],
                )
            )
            impactx_beamline.append(elements.Drift(name=name, ds=ds, nslice=nslice))
            impactx_beamline.append(
                elements.Aperture(
                    name=name + "_aperture_exit",
                    aperture_x=aperture_params["aperture_x"],
                    aperture_y=aperture_params["aperture_y"],
                    shape=aperture_params["shape"],
                    dx=aperture_params["dx"],
                    dy=aperture_params["dy"],
                )
            )
            return "aperture_drift_aperture"
        if ds > 0.0:
            impactx_beamline.append(elements.Drift(name=name, ds=ds, nslice=nslice))
            return "drift_only"
        if aperture_params is not None:
            impactx_beamline.append(elements.Aperture(**aperture_params))
            return "aperture_only"
        impactx_beamline.append(elements.Marker(name=name))
        return "marker_fallback"

    def make_bend_body_element(
        *,
        name,
        ds,
        rc,
        tilt_degree,
        k1,
        k1s,
        k2,
        k2s,
        k3,
        k3s,
    ):
        """Create the most faithful currently available ImpactX bend body model.

        Priority:
        1) ExactCFbend when skew or higher body multipoles are present
        2) CFbend for pure dipole+quadrupole
        3) Sbend for pure dipole geometry
        """
        use_exact_cfbend = any(abs(v) > 0.0 for v in (k1s, k2, k2s, k3, k3s))
        if use_exact_cfbend:
            curvature = 1.0 / rc if rc != 0.0 else 0.0
            return elements.ExactCFbend(
                name=name,
                ds=ds,
                k_normal=[curvature, k1, k2, k3],
                k_skew=[0.0, k1s, k2s, k3s],
                unit=0,
                rotation=tilt_degree,
                nslice=nslice,
            )
        if abs(k1) > 0.0:
            return elements.CFbend(
                name=name,
                ds=ds,
                rc=rc,
                k=k1,
                rotation=tilt_degree,
                nslice=nslice,
            )
        return elements.Sbend(
            name=name,
            ds=ds,
            rc=rc,
            rotation=tilt_degree,
            nslice=nslice,
        )

    for d_raw in parsed_beamline:
        d = _TrackedElement(d_raw)
        # print(d)
        if d["type"] in [k.casefold() for k in list(madx_to_impactx_dict.keys())]:
            if d["type"] == "marker":
                # Always preserve MARKER semantics as an explicit no-op element.
                impactx_beamline.append(elements.Marker(name=d["name"]))

                # MARKER can additionally carry aperture metadata in MAD-X.
                # Map supported aperture descriptions to an explicit ImpactX Aperture element.
                aperture_element = marker_aperture_element(d)
                if aperture_element is not None:
                    impactx_beamline.append(aperture_element)
                elif d.get("apertype", "") or d.get("aperture", []):
                    # TODO(audit): Add mapping for additional MAD-X aperture types
                    # (e.g. rectellipse/racetrack/octagon and vertex-defined profiles).
                    _warn(
                        f"MARKER '{d['name']}' aperture metadata could not be mapped exactly (apertype='{d.get('apertype', '')}'); skipping aperture translation."
                    )
            elif d["type"] in ("collimator", "rcollimator", "ecollimator"):
                ds = d.get("l", 0.0)
                aperture_params = aperture_params_from_madx(d)
                # MAD-X legacy XSIZE/YSIZE are documented obsolete and not used.
                # Read them explicitly so unread-attribute warnings stay meaningful.
                _ = d.get("xsize", 0.0)
                _ = d.get("ysize", 0.0)
                mode = append_collimator_equivalent(
                    name=d["name"], ds=ds, aperture_params=aperture_params
                )

                # TODO(audit): MAD-X collimator material/jaw interactions are not
                # represented in ImpactX; we preserve only optics and aperture checks.
                if mode == "aperture_drift_aperture":
                    _warn(
                        f"{d['type'].upper()} '{d['name']}' is approximated as aperture(entry)+drift+aperture(exit); jaw material interactions are not translated."
                    )
                elif mode == "drift_only":
                    _warn(
                        f"{d['type'].upper()} '{d['name']}' has no translatable aperture geometry; using Drift only."
                    )
                elif mode == "aperture_only":
                    _warn(
                        f"{d['type'].upper()} '{d['name']}' is modeled as a zero-length Aperture; jaw material interactions are not translated."
                    )
                else:
                    _warn(
                        f"{d['type'].upper()} '{d['name']}' has neither finite length nor translatable aperture geometry; using Marker fallback."
                    )
            elif d["type"] == "placeholder":
                ds = d.get("l", d.get("lrad", 0.0))
                if ds > 0.0:
                    impactx_beamline.append(
                        elements.Drift(name=d["name"], ds=ds, nslice=nslice)
                    )
                    _warn(
                        f"PLACEHOLDER '{d['name']}' is approximated as Drift(L={ds})."
                    )
                else:
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    _warn(f"PLACEHOLDER '{d['name']}' has zero length; using Marker.")
            elif d["type"] == "instrument":
                # MAD-X manual 11.20: INSTRUMENT is optically drift-like.
                ds = d.get("l", 0.0)
                if ds > 0.0:
                    impactx_beamline.append(
                        elements.Drift(name=d["name"], ds=ds, nslice=nslice)
                    )
                else:
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    _warn(f"INSTRUMENT '{d['name']}' has zero length; using Marker.")
            elif d["type"] == "drift":
                ds = d.get("l", 0.0)
                if ds > 0:
                    impactx_beamline.append(
                        elements.Drift(name=d["name"], ds=ds, nslice=nslice)
                    )
            elif d["type"] == "quadrupole":
                ds = d.get("l", 0.0)
                k1 = d.get("k1", 0.0)
                k1s = d.get("k1s", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if ds > 0:
                    if abs(k1s) > 0:
                        impactx_beamline.append(
                            elements.ExactMultipole(
                                name=d["name"],
                                ds=ds,
                                k_normal=[0.0, k1],
                                k_skew=[0.0, k1s],
                                rotation=tilt_degree,
                                nslice=nslice,
                            )
                        )
                    else:
                        impactx_beamline.append(
                            elements.Quad(
                                name=d["name"],
                                ds=ds,
                                k=k1,
                                rotation=tilt_degree,
                                nslice=nslice,
                            )
                        )
                else:
                    # Thin quad: integrated strength KL = k1 * l
                    knl = d.get("knl", [])
                    ksl = d.get("ksl", [])
                    k1l = knl[1] if len(knl) > 1 else k1 * ds
                    k1sl = ksl[1] if len(ksl) > 1 else k1s * ds
                    impactx_beamline.append(
                        elements.Multipole(
                            name=d["name"],
                            # MAD-X multipole order n maps to ImpactX multipole index n+1:
                            # n=1 (quadrupole) -> multipole=2.
                            multipole=2,
                            K_normal=k1l,
                            K_skew=k1sl,
                            rotation=tilt_degree,
                        )
                    )
            elif d["type"] == "sbend":
                ds = d.get("l", 0.0)
                angle = d.get("angle", 0.0)
                k1 = d.get("k1", 0.0)
                k1s = d.get("k1s", 0.0)
                k2 = d.get("k2", 0.0)
                k2s = d.get("k2s", 0.0)
                k3 = d.get("k3", 0.0)
                k3s = d.get("k3s", 0.0)
                e1 = d.get("e1", 0.0)
                e2 = d.get("e2", 0.0)
                hgap = d.get("hgap", 0.0)
                fint = d.get("fint", 0.0)
                fintx = d.get("fintx", fint)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if ds > 0:
                    if abs(angle) > 0:
                        rc = ds / angle
                        has_edges = (
                            abs(e1) > 0
                            or abs(e2) > 0
                            or abs(hgap) > 0
                            or abs(fint) > 0
                            or abs(fintx) > 0
                        )

                        if has_edges:
                            impactx_beamline.append(
                                elements.DipEdge(
                                    name=d["name"] + "_edge_entry",
                                    psi=e1,
                                    rc=rc,
                                    g=2.0 * hgap,
                                    K2=fint,
                                    location="entry",
                                    rotation=tilt_degree,
                                )
                            )

                        impactx_beamline.append(
                            make_bend_body_element(
                                name=d["name"],
                                ds=ds,
                                rc=rc,
                                tilt_degree=tilt_degree,
                                k1=k1,
                                k1s=k1s,
                                k2=k2,
                                k2s=k2s,
                                k3=k3,
                                k3s=k3s,
                            )
                        )

                        if has_edges:
                            impactx_beamline.append(
                                elements.DipEdge(
                                    name=d["name"] + "_edge_exit",
                                    psi=e2,
                                    rc=rc,
                                    g=2.0 * hgap,
                                    K2=fintx,
                                    location="exit",
                                    rotation=tilt_degree,
                                )
                            )
                    else:
                        # TODO(audit): MAD-X allows ANGLE=0 SBEND with extra attributes.
                        # ImpactX has no direct "straight CF bend" with all SBEND options.
                        # Fallback to Quad or Drift and warn.
                        if abs(k1) > 0:
                            impactx_beamline.append(
                                elements.Quad(
                                    name=d["name"],
                                    ds=ds,
                                    k=k1,
                                    rotation=tilt_degree,
                                    nslice=nslice,
                                )
                            )
                        else:
                            impactx_beamline.append(
                                elements.Drift(name=d["name"], ds=ds, nslice=nslice)
                            )
                        _warn(
                            f"SBEND '{d['name']}' has ANGLE=0; using simplified straight-element fallback."
                        )
                elif abs(angle) > 0:
                    # Thin dipole kick
                    if any(abs(v) > 0.0 for v in (k1, k1s, k2, k2s, k3, k3s, e1, e2)):
                        # TODO(audit): Thin SBEND with combined-function/edge effects is
                        # currently approximated as a pure thin dipole in ImpactX.
                        _warn(
                            f"Thin SBEND '{d['name']}' has unsupported non-dipole features; using ThinDipole fallback."
                        )
                    impactx_beamline.append(
                        elements.ThinDipole(name=d["name"], theta=angle, rc=0.0)
                    )
                for attr in unsupported_bend_attrs:
                    if abs(d.get(attr, 0.0)) > 0:
                        # TODO(audit): Map higher-order/extended SBEND attributes when
                        # ImpactX provides a direct equivalent in this translator path.
                        _warn(
                            f"SBEND '{d['name']}' attribute '{attr}' is not translated; using simplified model."
                        )
            elif d["type"] == "rbend":
                # RBEND: rectangular bending magnet
                # Converts to DipEdge + SBend + DipEdge
                # Reference: https://mad.web.cern.ch/mad/webguide/manual.html#Ch11.S3
                # TODO(audit): Confirm this decomposition preserves MAD-X RBEND
                # physics for all relevant options (RBARC, edge/fringe conventions,
                # and possible combined-function strengths).
                angle = d.get("angle", 0.0)
                l_chord = d.get("l", 0.0)
                k1 = d.get("k1", 0.0)
                k1s = d.get("k1s", 0.0)
                k2 = d.get("k2", 0.0)
                k2s = d.get("k2s", 0.0)
                k3 = d.get("k3", 0.0)
                k3s = d.get("k3s", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg

                if l_chord > 0:
                    # RBARC option controls how RBEND L is interpreted:
                    # - rbarc=true  (default): L is straight/chord length -> convert to arc
                    # - rbarc=false: L is already arc length -> do not convert
                    if abs(angle) > 0 and rbarc:
                        l_arc = l_chord * angle / (2.0 * math.sin(angle / 2.0))
                    else:
                        l_arc = l_chord

                    rc = l_arc / angle if abs(angle) > 0 else 0.0

                    # RBEND edge angles: effective = E1/E2 + ANGLE/2
                    e1 = d.get("e1", 0.0) + angle / 2.0
                    e2 = d.get("e2", 0.0) + angle / 2.0

                    hgap = d.get("hgap", 0.0)
                    fint = d.get("fint", 0.0)
                    fintx = d.get("fintx", fint)

                    # Entry DipEdge
                    impactx_beamline.append(
                        elements.DipEdge(
                            name=d["name"] + "_edge_entry",
                            psi=e1,
                            rc=rc,
                            g=2.0 * hgap,
                            K2=fint,
                            location="entry",
                            rotation=tilt_degree,
                        )
                    )

                    # Body: prefer combined-function model if representable.
                    impactx_beamline.append(
                        make_bend_body_element(
                            name=d["name"],
                            ds=l_arc,
                            rc=rc,
                            tilt_degree=tilt_degree,
                            k1=k1,
                            k1s=k1s,
                            k2=k2,
                            k2s=k2s,
                            k3=k3,
                            k3s=k3s,
                        )
                    )
                    # Exit DipEdge
                    # Exit DipEdge
                    impactx_beamline.append(
                        elements.DipEdge(
                            name=d["name"] + "_edge_exit",
                            psi=e2,
                            rc=rc,
                            g=2.0 * hgap,
                            K2=fintx,
                            location="exit",
                            rotation=tilt_degree,
                        )
                    )
                elif abs(angle) > 0:
                    if (
                        any(abs(v) > 0.0 for v in (k1, k1s, k2, k2s, k3, k3s))
                        or abs(d.get("e1", 0.0)) > 0
                        or abs(d.get("e2", 0.0)) > 0
                    ):
                        # TODO(audit): Thin RBEND with combined-function/edge effects is
                        # currently approximated as a pure thin dipole in ImpactX.
                        _warn(
                            f"Thin RBEND '{d['name']}' has unsupported non-dipole features; using ThinDipole fallback."
                        )
                    impactx_beamline.append(
                        elements.ThinDipole(name=d["name"], theta=angle, rc=0.0)
                    )
                for attr in unsupported_bend_attrs:
                    if abs(d.get(attr, 0.0)) > 0:
                        # TODO(audit): Map higher-order/extended RBEND attributes when
                        # ImpactX provides a direct equivalent in this translator path.
                        _warn(
                            f"RBEND '{d['name']}' attribute '{attr}' is not translated; using simplified model."
                        )
            elif d["type"] == "solenoid":
                ds = d.get("l", 0.0)
                ks = d.get("ks", 0.0)
                ksi = d.get("ksi", 0.0)
                lrad = d.get("lrad", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if ds > 0:
                    impactx_beamline.append(
                        elements.Sol(
                            name=d["name"],
                            ds=ds,
                            ks=ks,
                            rotation=tilt_degree,
                            nslice=nslice,
                        )
                    )
                    if abs(ksi) > 0:
                        # TODO(audit): Clarify precedence/combination rules for KS and KSI
                        # in thick MAD-X SOLENOID beyond current direct KS mapping.
                        _warn(
                            f"SOLENOID '{d['name']}' provides KSI with L>0; using KS mapping only."
                        )
                elif abs(ksi) > 0 and lrad > 0 and abs(ks) == 0:
                    # TODO(audit): MAD-X supports thin solenoid via KSI and LRAD.
                    # ImpactX has no dedicated thin solenoid element in this path.
                    # Approximate with an equivalent-thickness Sol over LRAD.
                    _warn(
                        f"SOLENOID '{d['name']}' thin model (L=0, KSI, LRAD) approximated as thick Sol over LRAD."
                    )
                    impactx_beamline.append(
                        elements.Sol(
                            name=d["name"],
                            ds=lrad,
                            ks=ksi / lrad,
                            rotation=tilt_degree,
                            nslice=nslice,
                        )
                    )
                elif abs(ks) > 0 or abs(ksi) > 0:
                    # TODO(audit): Add a dedicated thin-solenoid translation once
                    # ImpactX representation is available in this converter.
                    # Corner case: MAD-X thin-solenoid configurations without a usable LRAD
                    # cannot preserve integrated rotation in the current ImpactX element set.
                    _warn(
                        f"SOLENOID '{d['name']}' has nonzero strength but no translatable finite-length model; skipping element."
                    )
            elif d["type"] == "dipedge":
                h = d.get("h", 0.0)
                # MAD-X stores an ENTRANCE logical on DIPEDGE (e.g., from MAKETHIN).
                # Map directly to ImpactX DipEdge location.
                location = "entry" if bool(d.get("entrance", False)) else "exit"
                impactx_beamline.append(
                    elements.DipEdge(
                        name=d["name"],
                        psi=d.get("e1", 0.0),
                        rc=1.0 / h if h != 0.0 else 0.0,
                        # MAD-X is using half the gap height
                        g=2.0 * d.get("hgap", 0.0),
                        K2=d.get("fint", 0.0),
                        location=location,
                        rotation=d.get("tilt", 0.0) * rad_to_deg,
                    )
                )
            elif d["type"] in ("kicker", "tkicker"):
                # Corner case: ImpactX Kicker is thin. For MAD-X finite L, we preserve
                # placement/length with a drift-kick-drift composition.
                ds = d.get("l", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if abs(d.get("sinkick", 0.0)) > 0:
                    # TODO(audit): Map sinusoidally driven kicker modes (SINKICK/SINPEAK/SINTUNE/SINPHASE).
                    _warn(
                        f"KICKER '{d['name']}' uses sinusoidal kick options not translated; using static kick values."
                    )
                append_thin_with_length(
                    elements.Kicker(
                        name=d["name"],
                        xkick=d.get("hkick", 0.0),
                        ykick=d.get("vkick", 0.0),
                        rotation=tilt_degree,
                    ),
                    ds,
                    d["name"],
                )
            elif d["type"] == "hkicker":
                ds = d.get("l", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if abs(d.get("sinkick", 0.0)) > 0:
                    # TODO(audit): Map sinusoidally driven kicker modes (SINKICK/SINPEAK/SINTUNE/SINPHASE).
                    _warn(
                        f"HKICKER '{d['name']}' uses sinusoidal kick options not translated; using static kick values."
                    )
                append_thin_with_length(
                    elements.Kicker(
                        name=d["name"],
                        xkick=d.get("kick", d.get("hkick", 0.0)),
                        ykick=0.0,
                        rotation=tilt_degree,
                    ),
                    ds,
                    d["name"],
                )
            elif d["type"] == "vkicker":
                ds = d.get("l", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if abs(d.get("sinkick", 0.0)) > 0:
                    # TODO(audit): Map sinusoidally driven kicker modes (SINKICK/SINPEAK/SINTUNE/SINPHASE).
                    _warn(
                        f"VKICKER '{d['name']}' uses sinusoidal kick options not translated; using static kick values."
                    )
                append_thin_with_length(
                    elements.Kicker(
                        name=d["name"],
                        xkick=0.0,
                        ykick=d.get("kick", d.get("vkick", 0.0)),
                        rotation=tilt_degree,
                    ),
                    ds,
                    d["name"],
                )
            elif d["type"] in ("monitor", "hmonitor", "vmonitor"):
                if d.get("l", 0.0) > 0:
                    impactx_beamline.append(
                        elements.Drift(
                            name=d["name"] + "_drift", ds=d["l"], nslice=nslice
                        )
                    )
                impactx_beamline.append(
                    elements.BeamMonitor(name=d["name"], backend="h5")
                )
            elif d["type"] == "sextupole":
                ds = d.get("l", 0.0)
                k2 = d.get("k2", 0.0)
                k2s = d.get("k2s", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if ds > 0:
                    impactx_beamline.append(
                        elements.ExactMultipole(
                            name=d["name"],
                            ds=ds,
                            k_normal=[0.0, 0.0, k2],
                            k_skew=[0.0, 0.0, k2s],
                            rotation=tilt_degree,
                            nslice=nslice,
                        )
                    )
                else:
                    # Thin branch: prefer integrated strengths from KNL/KSL arrays.
                    knl = d.get("knl", [])
                    ksl = d.get("ksl", [])
                    k2l = knl[2] if len(knl) > 2 else k2 * ds
                    k2sl = ksl[2] if len(ksl) > 2 else k2s * ds
                    impactx_beamline.append(
                        elements.Multipole(
                            name=d["name"],
                            # MAD-X n=2 (sextupole) -> ImpactX multipole=3.
                            multipole=3,
                            K_normal=k2l,
                            K_skew=k2sl,
                            rotation=tilt_degree,
                        )
                    )
            elif d["type"] == "octupole":
                ds = d.get("l", 0.0)
                k3 = d.get("k3", 0.0)
                k3s = d.get("k3s", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if ds > 0:
                    impactx_beamline.append(
                        elements.ExactMultipole(
                            name=d["name"],
                            ds=ds,
                            k_normal=[0.0, 0.0, 0.0, k3],
                            k_skew=[0.0, 0.0, 0.0, k3s],
                            rotation=tilt_degree,
                            nslice=nslice,
                        )
                    )
                else:
                    # Thin branch: prefer integrated strengths from KNL/KSL arrays.
                    knl = d.get("knl", [])
                    ksl = d.get("ksl", [])
                    k3l = knl[3] if len(knl) > 3 else k3 * ds
                    k3sl = ksl[3] if len(ksl) > 3 else k3s * ds
                    impactx_beamline.append(
                        elements.Multipole(
                            name=d["name"],
                            # MAD-X n=3 (octupole) -> ImpactX multipole=4.
                            multipole=4,
                            K_normal=k3l,
                            K_skew=k3sl,
                            rotation=tilt_degree,
                        )
                    )
            elif d["type"] == "multipole":
                # MAD-X MULTIPOLE is a thin element with KNL/KSL arrays
                # KNL = integrated normal strengths, KSL = integrated skew strengths
                knl = d.get("knl", [])
                ksl = d.get("ksl", [])
                ds = d.get("l", 0.0)
                lrad = d.get("lrad", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                max_order = max(len(knl), len(ksl)) - 1

                if ds > 0.0:
                    # TODO(audit): MAD-X MULTIPOLE is conceptually thin. Finite length is often
                    # used for weak-focusing/radiation modeling (e.g., via LRAD). We preserve
                    # length with drift-kick-drift and warn.
                    _warn(
                        f"MULTIPOLE '{d['name']}' has finite L={ds}; using drift-kick-drift approximation."
                    )

                if lrad > 0.0:
                    # TODO(audit): LRAD/thin_foc behavior has no direct translation path here.
                    if thin_foc:
                        _warn(
                            f"MULTIPOLE '{d['name']}' uses LRAD={lrad} with OPTION,THIN_FOC=true; weak-focusing/radiation modeling is not translated."
                        )
                    else:
                        _warn(
                            f"MULTIPOLE '{d['name']}' uses LRAD={lrad}; LRAD-related modeling is not translated."
                        )

                if max_order >= 0:
                    if ds > 0.0:
                        impactx_beamline.append(
                            elements.Drift(
                                name=d["name"] + "_drift_in",
                                ds=0.5 * ds,
                                nslice=nslice,
                            )
                        )
                    for order in range(max_order + 1):
                        kn = knl[order] if order < len(knl) else 0.0
                        ks = ksl[order] if order < len(ksl) else 0.0
                        if kn != 0.0 or ks != 0.0:
                            impactx_beamline.append(
                                elements.Multipole(
                                    name=d["name"],
                                    # MAD-X MULTIPOLE arrays are 0-based in order:
                                    # index 0=dipole, 1=quadrupole, ... ; ImpactX uses
                                    # 1=dipole, 2=quadrupole, ... hence (order + 1).
                                    multipole=order + 1,
                                    K_normal=kn,
                                    K_skew=ks,
                                    rotation=tilt_degree,
                                )
                            )
                    if ds > 0.0:
                        impactx_beamline.append(
                            elements.Drift(
                                name=d["name"] + "_drift_out",
                                ds=0.5 * ds,
                                nslice=nslice,
                            )
                        )
            elif d["type"] == "nllens":
                if abs(d.get("cnll", 0.0)) == 0.0:
                    # TODO(audit): CNLL=0 is singular in MAD-X NLLENS model.
                    _warn(
                        f"NLLENS '{d['name']}' has CNLL=0; model may be singular/undefined."
                    )
                impactx_beamline.append(
                    elements.NonlinearLens(
                        name=d["name"],
                        knll=d.get("knll", 0.0),
                        cnll=d.get("cnll", 0.0),
                    )
                )
            elif d["type"] == "rfcavity":
                # MAD-X: volt [MV], lag [2pi], harmon [1]
                # ImpactX ShortRF: V [MV], freq [Hz], phase [deg]
                # ImpactX RFCavity: escale [MV/m], freq [Hz], phase [deg]
                # TODO(audit): Verify sign/unit conventions for LAG->phase conversion
                # and dependency on BEAM FREQ0 in MAD-X implementation.
                volt = d.get("volt", 0.0)  # MV
                lag = d.get("lag", 0.0)  # fraction of 2pi
                harmon = d.get("harmon", 1.0)
                phase = lag * 360.0  # convert to degrees
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if "freq" in d and d.get("freq", 0.0) > 0.0:
                    # MAD-X RFCAVITY frequency is specified in MHz.
                    freq = d.get("freq", 0.0) * 1.0e6
                else:
                    freq = harmon * freq0 * 1.0e6  # MHz -> Hz
                    if freq0 == 0.0:
                        # Corner case: neither FREQ nor BEAM FREQ0 are usable.
                        # We keep translation deterministic (freq=0) and warn.
                        _warn(
                            f"RFCAVITY '{d['name']}' has no explicit FREQ and BEAM FREQ0=0; resulting RF frequency is 0."
                        )
                ds = d.get("l", 0.0)

                for attr in (
                    "betrf",
                    "pg",
                    "shunt",
                    "tfill",
                    "n_bessel",
                    "no_cavity_totalpath",
                ):
                    if attr in d and abs(d.get(attr, 0.0)) > 0:
                        # TODO(audit): Extend translation for advanced MAD-X/PTC cavity attributes.
                        _warn(
                            f"RFCAVITY '{d['name']}' attribute '{attr}' is not translated; using simplified model."
                        )

                if ds > 0:
                    # Thick cavity: use RFCavity with pillbox field profile
                    escale = volt / ds  # MV/m
                    impactx_beamline.append(
                        elements.RFCavity(
                            name=d["name"],
                            ds=ds,
                            escale=escale,
                            freq=freq,
                            phase=phase,
                            cos_coefficients=[1.0],
                            sin_coefficients=[0.0],
                            rotation=tilt_degree,
                            nslice=nslice,
                        )
                    )
                else:
                    # Thin cavity: use ShortRF
                    impactx_beamline.append(
                        elements.ShortRF(
                            name=d["name"],
                            V=volt,
                            freq=freq,
                            phase=phase,
                            rotation=tilt_degree,
                        )
                    )
            warn_unread_element_attributes(d)
        else:
            raise NotImplementedError(
                "The beamline element named ",
                d["name"],
                "of type ",
                d["type"],
                "is not implemented in impactx.elements.",
                "Available elements are:",
                list(madx_to_impactx_dict.keys()),
            )
    return impactx_beamline


def read_lattice(madx_file, nslice=1):
    """
    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param madx_file: file name to MAD-X file with beamline elements
    :param nslice: number of ds slices per element
    :return: list of ImpactX.KnownElements
    """
    madx = MADXParser()
    madx.parse(madx_file)
    beamline = madx.getBeamline()
    freq0 = madx.getFreq0()
    all_options = madx.getOptions()
    options = {
        # MAD-X options are numeric in parser context; bool() captures 0/1 and expressions.
        "rbarc": bool(madx.getOption("rbarc", 1.0)),
        "thin_foc": bool(madx.getOption("thin_foc", 1.0)),
    }
    for opt_name in sorted(all_options.keys()):
        if opt_name not in options:
            _warn(
                f"OPTION '{opt_name}' was parsed but is not consumed by MAD-X->ImpactX translation."
            )

    return lattice(beamline, nslice, freq0=freq0, options=options)


def beam(particle, charge=None, mass=None, energy=None):
    """
    Function that converts a list of beam parameter dictionaries in the MADXParser format into ImpactX format

    Rules following https://mad.web.cern.ch/mad/releases/5.02.08/madxuguide.pdf pages 55f.

    :param str particle: reference particle name
    :param float charge: particle charge (proton charge units)
    :param float mass: particle mass (electron masses)
    :param float energy: total particle energy (GeV)
        - MAD-X default: 1 GeV
    :return dict: dictionary containing particle and beam attributes in ImpactX units
    """

    GeV2MeV = 1.0e3
    kg2MeV = sc.c**2 / sc.electron_volt * 1.0e-6
    muon_mass = sc.physical_constants["electron-muon mass ratio"][0] / sc.m_e
    if energy is None:
        energy_MeV = 1.0e3  # MAD-X default is 1 GeV total particle energy
    else:
        energy_MeV = energy * GeV2MeV

    impactx_beam = {
        "positron": {"mass": sc.m_e * kg2MeV, "charge": 1.0},
        "electron": {"mass": sc.m_e * kg2MeV, "charge": -1.0},
        "proton": {"mass": sc.m_p * kg2MeV, "charge": 1.0},
        "antiproton": {"mass": sc.m_p * kg2MeV, "charge": -1.0},
        "posmuon": {
            "mass": muon_mass * kg2MeV,
            "charge": 1.0,
        },  # positively charged muon (anti-muon)
        "negmuon": {
            "mass": muon_mass * kg2MeV,
            "charge": -1.0,
        },  # negatively charged muon
        "ion": {"mass": sc.m_u * kg2MeV, "charge": 1.0},
        "generic": {"mass": mass, "charge": charge},
    }

    if particle not in impactx_beam.keys():
        _warn(f"Particle species name '{particle}' not in '{impactx_beam.keys()}'")
        print(
            "Choosing generic particle species, using provided `charge`, `mass` and `energy`."
        )
        _particle = "generic"
    else:
        _particle = particle
        print(
            f"Choosing provided particle species '{particle}', ignoring potentially provided `charge` and `mass` and setting defaults.",
        )

    reference_particle = impactx_beam[_particle]
    reference_particle["energy"] = energy_MeV

    return reference_particle


def read_beam(ref: RefPart, madx_file):
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

    warnings.warn(
        "Our MAD-X parser is under active development and provided as a preview. "
        "Please check any loaded MAD-X beams very carefully. Please report your "
        "experience and bugs on our issue tracker: "
        "https://github.com/BLAST-ImpactX/impactx/issues",
        RuntimeWarning,
        stacklevel=2,
    )

    madx = MADXParser()
    madx.parse(madx_file)

    ref_particle_dict = beam(
        madx.getParticle(),  # if particle species is known, mass, charge, and potentially energy are set to default
        # TODO MADX parser needs to extract charge if it's given,
        # TODO MADX parser needs to extract mass if it's given,
        energy=float(madx.getEtot()),  # MADX default energy is in GeV
    )

    ref.set_charge_qe(ref_particle_dict["charge"])
    ref.set_mass_MeV(ref_particle_dict["mass"])
    ref.set_kin_energy_MeV(ref_particle_dict["energy"] - ref_particle_dict["mass"])

    return ref
