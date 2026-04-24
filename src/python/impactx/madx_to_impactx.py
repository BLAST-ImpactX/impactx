#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import math
import warnings

from impactx import Map6x6, RefPart, elements

from .MADXParser import MADXParser

# Single-shape mapping from MAD-X APERTYPE to ImpactX Aperture.shape.
# Edit this table first when adding support for new MAD-X aperture types.
# Composite shapes (rectellipse, lhcscreen/rectcircle) are not in this table;
# they are expanded into two back-to-back ImpactX Apertures by
# aperture_params_from_madx (intersection of a rectangle and an ellipse/circle).
APERTURE_SHAPE_MAP = {
    "rectangle": "rectangular",
    "rectangular": "rectangular",
    "ellipse": "elliptical",
    "circle": "elliptical",
}


class MADXImpactXTranslatorWarning(UserWarning):
    """Warning category for MAD-X to ImpactX translation fallbacks."""

    pass


def _warn(message: str):
    """Emit a translation warning without mutating process-wide warning formatting."""
    warnings.warn(message, MADXImpactXTranslatorWarning, stacklevel=2)


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


def lattice(
    parsed_beamline,
    nslice=1,
    freq0=0.0,
    options=None,
    ref_mass_MeV=None,
    bv=1.0,
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

    supported_madx_elements = {
        "MARKER",
        "BEAMBEAM",  # beam-beam interaction, no ImpactX equivalent
        "DRIFT",
        "SBEND",  # Sector Bending Magnet
        "RBEND",  # Rectangular Bending Magnet -> DipEdge + Sbend + DipEdge
        "SOLENOID",  # Ideal, thick Solenoid: MAD-X user guide 10.9 p78
        "QUADRUPOLE",
        "DIPEDGE",
        # Kicker, idealized thin element,
        # MADX also defines length "L" and a roll angle around the longitudinal axis "TILT"
        # https://mad.web.cern.ch/mad/webguide/manual.html#Ch11.S11
        "KICKER",
        "TKICKER",  # thin kicker, treated same as KICKER
        "HKICKER",  # horizontal kicker
        "VKICKER",  # vertical kicker
        # note: in MAD-X, this keeps track only of the beam centroid,
        # "In addition it serves to record the beam position for closed orbit correction."
        "MONITOR",  # drift + output diagnostics
        "HMONITOR",  # horizontal monitor
        "VMONITOR",  # vertical monitor
        "MULTIPOLE",
        "SEXTUPOLE",
        "OCTUPOLE",
        "RFCAVITY",
        "NLLENS",
        "COLLIMATOR",  # fallback: geometry only (+ drift for length)
        "RCOLLIMATOR",  # rectangular collimator alias
        "ECOLLIMATOR",  # elliptical collimator alias
        "PLACEHOLDER",  # non-physical placement/device placeholder
        "INSTRUMENT",  # optically drift-like diagnostics placeholder
        # TODO Figure out how to identify these
        "ShortRF",
        "ConstF",
    }
    supported_madx_element_types = {name.casefold() for name in supported_madx_elements}

    impactx_beamline = []
    options = options or {}
    # MAD-X defaults from OPTION command table:
    # rbarc=true, thin_foc=true
    rbarc = bool(options.get("rbarc", True))
    # TODO(audit): OPTION,THIN_FOC toggles MAD-X's weak-focusing contribution on
    # thin elements with LRAD>0 (mainly MULTIPOLE). We currently surface the flag
    # in MULTIPOLE warnings but do not translate a weak-focusing term into the
    # output lattice; doing so would require composing a per-MULTIPOLE thin
    # quadrupole kick of strength ~angle^2/LRAD when THIN_FOC=true. Track when
    # ImpactX grows a direct equivalent.
    thin_foc = bool(options.get("thin_foc", True))
    rad_to_deg = 180.0 / math.pi
    # MAD-X bend maps are driven by ANGLE/TILT; K0/K0S are obsolete for map construction
    # and should be consumed silently if present. We read them via
    # ``_ = d.get(attr, 0.0)`` in the SBEND/RBEND branches purely to mark them
    # "accessed" in the _TrackedElement wrapper, which suppresses the
    # unread-attribute warning that would otherwise fire.
    ignored_bend_attrs = ("k0", "k0s")
    # H1/H2 (pole-face curvature) still affect MAD-X fringe modeling and are not translated here.
    unsupported_bend_attrs = ("h1", "h2")
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
            # Note: this emission uses drift(L/2) + kick + drift(L/2), while MAD-X
            # uses the reverse symmetric split half-kick + drift(L) + half-kick
            # (twiss.f90:4587-4617). Both are second-order symmetric compositions
            # of drift + kick and are linearly identical; they differ only at
            # O(kick^2) (nonlinear) terms.
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
        """Build list of Aperture constructor kwargs from MAD-X aperture metadata.

        Returns a list (possibly empty) of kwargs dicts without a ``name`` field.
        Callers assign a context-appropriate name suffix.

        Simple MAD-X shapes (``rectangle``, ``ellipse``, ``circle``) yield a
        single-entry list. Composite shapes are modelled *exactly* as the
        intersection of two ImpactX Aperture primitives in series:
        - ``rectellipse`` -> [rectangular, elliptical]
        - ``lhcscreen`` / ``rectcircle`` -> [rectangular, elliptical-as-circle]
        (mad_aper.c:470-599).

        Returns:
            list of dicts with Aperture constructor kwargs (no ``name``). Empty
            list when the element has no translatable aperture metadata.
        """
        apertype = str(elem.get("apertype", "")).lower()
        aperture = elem.get("aperture", [])
        aper_offset = elem.get("aper_offset", [0.0, 0.0])

        if not (apertype or aperture):
            return []

        dx, dy = parse_aperture_offsets(aper_offset)
        # Aperture orientation in the lab frame is the sum of the element tilt
        # (which rotates the element's local frame) and aper_tilt (an additional
        # rotation of the aperture within that frame). Both default to 0 in MAD-X.
        rotation = (elem.get("tilt", 0.0) + elem.get("aper_tilt", 0.0)) * rad_to_deg
        base = {"dx": dx, "dy": dy, "rotation": rotation}

        # Composite shapes (intersection of rectangle and ellipse/circle).
        if apertype in ("rectellipse", "lhcscreen", "rectcircle"):
            if not isinstance(aperture, (list, tuple)) or len(aperture) < 3:
                return []
            rect_x = float(aperture[0])
            rect_y = float(aperture[1])
            if apertype == "rectellipse":
                if len(aperture) < 4:
                    return []
                ell_x = float(aperture[2])
                ell_y = float(aperture[3])
            else:
                # lhcscreen / rectcircle: aperture[2] is the circle radius
                ell_x = float(aperture[2])
                ell_y = float(aperture[2])
            if rect_x <= 0.0 or rect_y <= 0.0 or ell_x <= 0.0 or ell_y <= 0.0:
                return []
            return [
                {
                    "aperture_x": rect_x,
                    "aperture_y": rect_y,
                    "shape": "rectangular",
                    **base,
                },
                {
                    "aperture_x": ell_x,
                    "aperture_y": ell_y,
                    "shape": "elliptical",
                    **base,
                },
            ]

        # Single-shape mapping.
        shape = APERTURE_SHAPE_MAP.get(apertype)
        ax, ay = parse_aperture_half_axes(aperture)

        # Circle in MAD-X uses one radius; mirror to y if only x is provided.
        if apertype == "circle" and ay == 0.0:
            ay = ax

        if shape is None or ax <= 0.0 or ay <= 0.0:
            return []

        return [
            {"aperture_x": ax, "aperture_y": ay, "shape": shape, **base},
        ]

    def _aperture_name(base_name, suffix, params, total):
        """Name an Aperture element; disambiguate shape for composite shapes."""
        if total <= 1:
            return base_name + suffix
        # rectangular -> _rect, elliptical -> _ell
        tag = "_rect" if params["shape"] == "rectangular" else "_ell"
        return base_name + suffix + tag

    def aperture_elements_from_madx(elem, *, name_suffix="_aperture"):
        """Build ImpactX Aperture elements from MAD-X aperture metadata."""
        params_list = aperture_params_from_madx(elem)
        total = len(params_list)
        return [
            elements.Aperture(
                name=_aperture_name(elem["name"], name_suffix, p, total), **p
            )
            for p in params_list
        ]

    def marker_aperture_elements(marker):
        """Build ImpactX Aperture elements from MARKER aperture metadata."""
        return aperture_elements_from_madx(marker)

    def append_collimator_equivalent(*, name, ds, aperture_params_list):
        """Append a MAD-X collimator equivalent using current ImpactX primitives.

        MAD-X semantics (manual 11.21 and tracking source):
        - Optically behaves like a drift.
        - Aperture checks are applied at entrance.

        ImpactX has no dedicated collimator jaw-material model here, so we
        preserve drift optics and approximate the entrance check with a
        leading Aperture element (or two, for composite apertypes).
        """

        def _append_apertures(suffix):
            total = len(aperture_params_list)
            for p in aperture_params_list:
                impactx_beamline.append(
                    elements.Aperture(name=_aperture_name(name, suffix, p, total), **p)
                )

        has_aperture = bool(aperture_params_list)
        if ds > 0.0 and has_aperture:
            _append_apertures("_aperture_entry")
            impactx_beamline.append(elements.Drift(name=name, ds=ds, nslice=nslice))
            return "aperture_drift"
        if ds > 0.0:
            impactx_beamline.append(elements.Drift(name=name, ds=ds, nslice=nslice))
            return "drift_only"
        if has_aperture:
            _append_apertures("_aperture")
            return "aperture_only"
        impactx_beamline.append(elements.Marker(name=name))
        return "marker_fallback"

    def append_placeholder_equivalent(*, name, ds, aperture_params_list):
        """Append a MAD-X PLACEHOLDER equivalent while preserving aperture metadata."""

        def _append_apertures(suffix):
            total = len(aperture_params_list)
            for p in aperture_params_list:
                impactx_beamline.append(
                    elements.Aperture(name=_aperture_name(name, suffix, p, total), **p)
                )

        has_aperture = bool(aperture_params_list)
        if ds > 0.0:
            if has_aperture:
                _append_apertures("_aperture_entry")
                impactx_beamline.append(elements.Drift(name=name, ds=ds, nslice=nslice))
                return "aperture_drift"
            impactx_beamline.append(elements.Drift(name=name, ds=ds, nslice=nslice))
            return "drift_only"

        impactx_beamline.append(elements.Marker(name=name))
        if has_aperture:
            _append_apertures("_aperture")
            return "marker_aperture"
        return "marker_only"

    def normalized_rf_voltage(volt_MV):
        """Convert MAD-X VOLT [MV] to ImpactX ShortRF normalized voltage."""
        if ref_mass_MeV is None or ref_mass_MeV <= 0.0:
            _warn(
                "RFCAVITY reference particle mass is unavailable; using unnormalized MAD-X VOLT value as a fallback."
            )
            return volt_MV
        return volt_MV / ref_mass_MeV

    def normalized_rf_escale(volt_MV, ds):
        """Convert MAD-X VOLT/L [MV/m] to ImpactX RFCavity escale [1/m]."""
        if ref_mass_MeV is None or ref_mass_MeV <= 0.0:
            _warn(
                "RFCAVITY reference particle mass is unavailable; using unnormalized MAD-X VOLT/L value as a fallback."
            )
            return volt_MV / ds
        return volt_MV / (ds * ref_mass_MeV)

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

    def make_thin_dipole_kick(*, name, angle, tilt_degree):
        """Create a thin dipole kick fallback for zero-length bends.

        ImpactX ThinDipole expects degrees and a finite curvature radius, so a
        thin Multipole dipole is the safer direct representation for a MAD-X
        zero-length bend kick given in radians.
        """
        return elements.Multipole(
            name=name,
            multipole=1,
            K_normal=angle,
            K_skew=0.0,
            rotation=tilt_degree,
        )

    def normalized_fintx(*, fint, fintx_raw):
        """Map MAD-X FINTX sentinel values to the effective exit fringe integral."""
        if fintx_raw is None or fintx_raw < 0.0:
            return fint
        return fintx_raw

    # Element types whose per-type branch below handles aperture metadata
    # directly (MARKER merges Aperture with its no-op emission; collimators and
    # PLACEHOLDER use dedicated helpers). Every other supported type reuses the
    # uniform pre-dispatch emission below.
    aperture_handled_inline = (
        "marker",
        "collimator",
        "rcollimator",
        "ecollimator",
        "placeholder",
    )

    def emit_aperture_entry(elem):
        """Pre-emit ImpactX Aperture element(s) for generic MAD-X aperture metadata.

        MAD-X checks aperture at element entrance (trrun.f90:820-821 + trcoll),
        so we emit the Aperture primitive(s) *before* the main element. Composite
        MAD-X apertypes (rectellipse/lhcscreen/rectcircle) expand to two
        back-to-back ImpactX Apertures (rectangular + elliptical).
        """
        params_list = aperture_params_from_madx(elem)
        if params_list:
            total = len(params_list)
            for p in params_list:
                impactx_beamline.append(
                    elements.Aperture(
                        name=_aperture_name(elem["name"], "_aperture_entry", p, total),
                        **p,
                    )
                )
        elif elem.get("apertype", "") or elem.get("aperture", []):
            # TODO(audit): Add mapping for remaining MAD-X aperture types
            # (e.g. racetrack/octagon and vertex-defined profiles).
            _warn(
                f"{elem.get('type', '').upper()} '{elem['name']}' aperture metadata "
                f"could not be mapped exactly (apertype='{elem.get('apertype', '')}'); "
                f"skipping aperture translation."
            )

    for d_raw in parsed_beamline:
        d = _TrackedElement(d_raw)
        # print(d)
        if d["type"] in supported_madx_element_types:
            if d["type"] not in aperture_handled_inline:
                emit_aperture_entry(d)
            if d["type"] == "marker":
                # Always preserve MARKER semantics as an explicit no-op element.
                impactx_beamline.append(elements.Marker(name=d["name"]))

                # MARKER can additionally carry aperture metadata in MAD-X.
                # Map supported aperture descriptions to explicit ImpactX Aperture
                # element(s); composite apertypes (rectellipse/lhcscreen) expand to two.
                aperture_elements = marker_aperture_elements(d)
                if aperture_elements:
                    impactx_beamline.extend(aperture_elements)
                elif d.get("apertype", "") or d.get("aperture", []):
                    # TODO(audit): Add mapping for remaining MAD-X aperture types
                    # (e.g. racetrack/octagon and vertex-defined profiles).
                    _warn(
                        f"MARKER '{d['name']}' aperture metadata could not be mapped exactly (apertype='{d.get('apertype', '')}'); skipping aperture translation."
                    )
            elif d["type"] in ("collimator", "rcollimator", "ecollimator"):
                ds = d.get("l", 0.0)
                aperture_params_list = aperture_params_from_madx(d)
                # MAD-X legacy XSIZE/YSIZE are documented obsolete and not used.
                # Read them explicitly so unread-attribute warnings stay meaningful.
                _ = d.get("xsize", 0.0)
                _ = d.get("ysize", 0.0)
                mode = append_collimator_equivalent(
                    name=d["name"], ds=ds, aperture_params_list=aperture_params_list
                )

                # TODO(audit): MAD-X collimator material/jaw interactions are not
                # represented in ImpactX; we preserve only optics and aperture checks.
                if mode == "aperture_drift":
                    _warn(
                        f"{d['type'].upper()} '{d['name']}' is approximated as aperture(entry)+drift; jaw material interactions are not translated."
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
                aperture_params_list = aperture_params_from_madx(d)
                mode = append_placeholder_equivalent(
                    name=d["name"], ds=ds, aperture_params_list=aperture_params_list
                )
                if mode == "aperture_drift":
                    _warn(
                        f"PLACEHOLDER '{d['name']}' is approximated as aperture(entry)+drift to preserve its aperture metadata."
                    )
                elif mode == "drift_only":
                    _warn(
                        f"PLACEHOLDER '{d['name']}' is approximated as Drift(L={ds})."
                    )
                elif mode == "marker_aperture":
                    _warn(
                        f"PLACEHOLDER '{d['name']}' has zero length; using Marker + Aperture to preserve its aperture metadata."
                    )
                else:
                    _warn(f"PLACEHOLDER '{d['name']}' has zero length; using Marker.")
            elif d["type"] == "beambeam":
                impactx_beamline.append(elements.Marker(name=d["name"]))
                _warn(
                    f"BEAMBEAM '{d['name']}' has no ImpactX equivalent; using Marker fallback."
                )
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
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    if any(abs(v) > 0.0 for v in (k1, k1s, tilt_degree)):
                        _warn(
                            f"QUADRUPOLE '{d['name']}' has L=0 and no integrated-strength form in MAD-X; using Marker."
                        )
            elif d["type"] == "sbend":
                ds = d.get("l", 0.0)
                # MAD-X tmbend (twiss.f90:3856-3902) multiplies bvk into angle and
                # the body multipoles (k1, k1s, k2). Mirror that here. K0/K0S are
                # field-error placeholders MAD-X treats as obsolete; we already
                # drop them, so no bvk needed there.
                angle = bv * d.get("angle", 0.0)
                k0 = d.get("k0", 0.0)
                k0s = d.get("k0s", 0.0)
                k1 = bv * d.get("k1", 0.0)
                k1s = bv * d.get("k1s", 0.0)
                k2 = bv * d.get("k2", 0.0)
                k2s = bv * d.get("k2s", 0.0)
                k3 = bv * d.get("k3", 0.0)
                k3s = bv * d.get("k3s", 0.0)
                e1 = d.get("e1", 0.0)
                e2 = d.get("e2", 0.0)
                hgap = d.get("hgap", 0.0)
                fint = d.get("fint", 0.0)
                fintx = normalized_fintx(fint=fint, fintx_raw=d.get("fintx", None))
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if ds > 0:
                    if abs(angle) > 0:
                        rc = ds / angle
                        has_edges = (
                            abs(e1) > 0
                            or abs(e2) > 0
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
                        # MAD-X allows ANGLE=0 SBEND with extra attributes. Pick the
                        # narrowest ImpactX element that represents the remaining
                        # straight-field content (skew or higher multipoles require
                        # ExactMultipole).
                        if abs(k0) > 0.0 or abs(k0s) > 0.0:
                            _warn(
                                f"SBEND '{d['name']}' has ANGLE=0 but nonzero K0/K0S; these are dropped (MAD-X treats K0/K0S as obsolete)."
                            )
                        if any(abs(v) > 0.0 for v in (k1s, k2, k2s, k3, k3s)):
                            impactx_beamline.append(
                                elements.ExactMultipole(
                                    name=d["name"],
                                    ds=ds,
                                    k_normal=[0.0, k1, k2, k3],
                                    k_skew=[0.0, k1s, k2s, k3s],
                                    rotation=tilt_degree,
                                    nslice=nslice,
                                )
                            )
                        elif abs(k1) > 0:
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
                            f"SBEND '{d['name']}' has ANGLE=0; using straight-element fallback."
                        )
                elif abs(angle) > 0:
                    # Thin dipole kick
                    if any(abs(v) > 0.0 for v in (k1, k1s, k2, k2s, k3, k3s, e1, e2)):
                        # TODO(audit): Thin SBEND with combined-function/edge effects is
                        # currently approximated as a pure thin dipole kick in ImpactX.
                        _warn(
                            f"Thin SBEND '{d['name']}' has unsupported non-dipole features; using thin dipole kick fallback."
                        )
                    impactx_beamline.append(
                        make_thin_dipole_kick(
                            name=d["name"], angle=angle, tilt_degree=tilt_degree
                        )
                    )
                else:
                    # L=0 and ANGLE=0: preserve element as a Marker rather than
                    # silently dropping it from the lattice.
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    if any(abs(v) > 0.0 for v in (k1, k1s, k2, k2s, k3, k3s, e1, e2)):
                        _warn(
                            f"SBEND '{d['name']}' has L=0 and ANGLE=0 but nonzero strengths/edges; using Marker (content not translated)."
                        )
                    else:
                        _warn(f"SBEND '{d['name']}' has L=0 and ANGLE=0; using Marker.")
                for attr in ignored_bend_attrs:
                    _ = d.get(attr, 0.0)
                for attr in unsupported_bend_attrs:
                    if abs(d.get(attr, 0.0)) > 0:
                        # TODO(audit): Map higher-order/extended SBEND attributes when
                        # ImpactX provides a direct equivalent in this translator path.
                        _warn(
                            f"SBEND '{d['name']}' attribute '{attr}' is not translated; using simplified model."
                        )
            elif d["type"] == "rbend":
                # RBEND: rectangular bending magnet (parallel pole faces).
                # Decomposed here as DipEdge(E1+ANGLE/2) + SBend(L_arc) + DipEdge(E2+ANGLE/2).
                # The ANGLE/2 offset comes from the rectangular geometry itself: MAD-X
                # stores E1/E2 as user-specified (mad_elem.c:804 returns them as-is) and
                # applies the ANGLE/2 rectangular-geometry offset inside its RBEND tracking
                # kernels, not at attribute-read time. Manual ch. 11.3: "E1 and E2 refer
                # to the extra edge angles in addition to the half-angle of the RBEND."
                # Reference: https://mad.web.cern.ch/mad/webguide/manual.html#Ch11.S3
                # TODO(audit): Confirm this decomposition preserves MAD-X RBEND
                # physics for all relevant options (RBARC, edge/fringe conventions,
                # and possible combined-function strengths).
                # See SBEND above for bvk rationale (twiss.f90:3856-3902).
                angle = bv * d.get("angle", 0.0)
                l_chord = d.get("l", 0.0)
                k1 = bv * d.get("k1", 0.0)
                k1s = bv * d.get("k1s", 0.0)
                k2 = bv * d.get("k2", 0.0)
                k2s = bv * d.get("k2s", 0.0)
                k3 = bv * d.get("k3", 0.0)
                k3s = bv * d.get("k3s", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg

                if l_chord > 0:
                    # RBARC option controls how RBEND L is interpreted:
                    # - rbarc=true  (default): L is straight/chord length -> convert to arc
                    # - rbarc=false: L is already arc length -> do not convert
                    if abs(angle) > 0 and rbarc:
                        l_arc = l_chord * angle / (2.0 * math.sin(angle / 2.0))
                    else:
                        l_arc = l_chord

                    if abs(angle) > 0:
                        rc = l_arc / angle

                        # RBEND edge angles: effective = E1/E2 + ANGLE/2
                        e1 = d.get("e1", 0.0) + angle / 2.0
                        e2 = d.get("e2", 0.0) + angle / 2.0

                        hgap = d.get("hgap", 0.0)
                        fint = d.get("fint", 0.0)
                        fintx = normalized_fintx(
                            fint=fint, fintx_raw=d.get("fintx", None)
                        )

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
                        # ANGLE=0 RBEND: pick the narrowest straight-multipole element
                        # that represents the residual content. Skew or higher
                        # multipoles require ExactMultipole.
                        if any(abs(v) > 0.0 for v in (k1s, k2, k2s, k3, k3s)):
                            impactx_beamline.append(
                                elements.ExactMultipole(
                                    name=d["name"],
                                    ds=l_arc,
                                    k_normal=[0.0, k1, k2, k3],
                                    k_skew=[0.0, k1s, k2s, k3s],
                                    rotation=tilt_degree,
                                    nslice=nslice,
                                )
                            )
                        elif abs(k1) > 0:
                            impactx_beamline.append(
                                elements.Quad(
                                    name=d["name"],
                                    ds=l_arc,
                                    k=k1,
                                    rotation=tilt_degree,
                                    nslice=nslice,
                                )
                            )
                        else:
                            impactx_beamline.append(
                                elements.Drift(name=d["name"], ds=l_arc, nslice=nslice)
                            )
                        _warn(
                            f"RBEND '{d['name']}' has ANGLE=0; using straight-element fallback."
                        )
                elif abs(angle) > 0:
                    if (
                        any(abs(v) > 0.0 for v in (k1, k1s, k2, k2s, k3, k3s))
                        or abs(d.get("e1", 0.0)) > 0
                        or abs(d.get("e2", 0.0)) > 0
                    ):
                        # TODO(audit): Thin RBEND with combined-function/edge effects is
                        # currently approximated as a pure thin dipole kick in ImpactX.
                        _warn(
                            f"Thin RBEND '{d['name']}' has unsupported non-dipole features; using thin dipole kick fallback."
                        )
                    impactx_beamline.append(
                        make_thin_dipole_kick(
                            name=d["name"], angle=angle, tilt_degree=tilt_degree
                        )
                    )
                else:
                    # L=0 and ANGLE=0: preserve element as a Marker rather than
                    # silently dropping it from the lattice.
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    if any(abs(v) > 0.0 for v in (k1, k1s, k2, k2s, k3, k3s)) or any(
                        abs(d.get(attr, 0.0)) > 0.0 for attr in ("e1", "e2")
                    ):
                        _warn(
                            f"RBEND '{d['name']}' has L=0 and ANGLE=0 but nonzero strengths/edges; using Marker (content not translated)."
                        )
                    else:
                        _warn(f"RBEND '{d['name']}' has L=0 and ANGLE=0; using Marker.")
                for attr in ignored_bend_attrs:
                    _ = d.get(attr, 0.0)
                for attr in unsupported_bend_attrs:
                    if abs(d.get(attr, 0.0)) > 0:
                        # TODO(audit): Map higher-order/extended RBEND attributes when
                        # ImpactX provides a direct equivalent in this translator path.
                        _warn(
                            f"RBEND '{d['name']}' attribute '{attr}' is not translated; using simplified model."
                        )
            elif d["type"] == "solenoid":
                ds = d.get("l", 0.0)
                # MAD-X tmsol/tmsol_th apply bvk to both ks and ksi
                # (twiss.f90:6852-6853 and 8534-8535).
                ks = bv * d.get("ks", 0.0)
                ksi = bv * d.get("ksi", 0.0)
                lrad = d.get("lrad", 0.0)
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if ds > 0:
                    if abs(lrad) > 0.0:
                        _warn(
                            f"SOLENOID '{d['name']}' attribute 'lrad' is not translated for thick solenoids; using KS/KSI-only model."
                        )
                    if abs(ks) == 0.0 and abs(ksi) > 0.0:
                        # MAD-X derives KS from KSI/L for thick solenoids.
                        ks = ksi / ds
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
                            f"SOLENOID '{d['name']}' provides both KS and KSI with L>0; using KS mapping."
                        )
                elif abs(ks) > 0 or abs(ksi) > 0:
                    # MAD-X defines a dedicated zero-length thin-solenoid map.
                    # Until ImpactX grows a matching thin-solenoid element, use a
                    # zero-length LinearMap instead of inventing a finite geometry.
                    # See ImpactX issue #1413.
                    #
                    # MAD-X TWISS tmsol_th first-order map at the reference particle:
                    #   sk  = KS/2, skl = KSI/2, Q0 = skl, Q = -sk*Q0
                    #   x_f  =  C x + S y
                    #   px_f = QC x + QS y + C px + S py
                    #   y_f  = -S x + C y
                    #   py_f = -QS x + QC y - S px + C py
                    #
                    # TODO(audit): MAD-X tmsol_th (twiss.f90:8544) uses
                    #   Q0 = skl / (1 + deltap)
                    # i.e. the effective rotation angle is momentum-dependent. A static
                    # 6x6 LinearMap cannot express per-particle (1+deltap) dependence
                    # for the transverse block. We therefore evaluate C, S, Q at
                    # deltap=0 (reference particle), which is exact on-momentum but
                    # neglects chromatic (off-momentum) defocusing of the thin
                    # solenoid. A proper fix requires a dedicated thin-solenoid
                    # element in ImpactX that evaluates per-particle (1+deltap);
                    # track with ImpactX issue #1413.
                    C = math.cos(0.5 * ksi)
                    S = math.sin(0.5 * ksi)
                    Q = -0.25 * ks * ksi

                    R = Map6x6.identity()
                    R[1, 1] = C
                    R[1, 2] = 0.0
                    R[1, 3] = S
                    R[1, 4] = 0.0
                    R[2, 1] = Q * C
                    R[2, 2] = C
                    R[2, 3] = Q * S
                    R[2, 4] = S
                    R[3, 1] = -S
                    R[3, 2] = 0.0
                    R[3, 3] = C
                    R[3, 4] = 0.0
                    R[4, 1] = -Q * S
                    R[4, 2] = -S
                    R[4, 3] = Q * C
                    R[4, 4] = C

                    impactx_beamline.append(
                        elements.LinearMap(
                            name=d["name"],
                            R=R,
                            ds=0.0,
                            rotation=tilt_degree,
                        )
                    )
                    _warn(
                        f"SOLENOID '{d['name']}' thin model (L=0) is approximated by a zero-length LinearMap from the MAD-X first-order map; nonlinear tracking terms and LRAD/radiation effects are not translated. See ImpactX issue #1413."
                    )
            elif d["type"] == "dipedge":
                h = d.get("h", 0.0)
                he = d.get("he", 0.0)
                # MAD-X stores an ENTRANCE logical on DIPEDGE (e.g., from MAKETHIN).
                # Map directly to ImpactX DipEdge location.
                # Note: MAD-X's ENTRANCE default is false (mad_dict.c:2957), so a
                # hand-written DIPEDGE without an explicit ENTRANCE=true is translated
                # as an exit edge -- this matches MAD-X's stored default but may be
                # surprising to users who treat standalone DIPEDGE as an entry edge.
                location = "entry" if bool(d.get("entrance", False)) else "exit"
                if abs(h) == 0.0:
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    if any(
                        abs(v) > 0.0
                        for v in (
                            d.get("e1", 0.0),
                            d.get("fint", 0.0),
                            d.get("hgap", 0.0),
                            he,
                            d.get("tilt", 0.0),
                        )
                    ):
                        _warn(
                            f"DIPEDGE '{d['name']}' has H=0 and is a no-op in MAD-X; using Marker."
                        )
                else:
                    if abs(he) > 0.0:
                        # MAD-X HE is the pole-face curvature. ImpactX DipEdge does not
                        # expose a direct HE input, so this higher-order contribution is dropped.
                        _warn(
                            f"DIPEDGE '{d['name']}' attribute 'he' is not translated; using simplified model."
                        )
                    impactx_beamline.append(
                        elements.DipEdge(
                            name=d["name"],
                            psi=d.get("e1", 0.0),
                            rc=1.0 / h,
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
                # bvk applies to the kicker kicks (twiss.f90:4561-4563).
                ds = d.get("l", 0.0)
                xkick = bv * (d.get("hkick", 0.0) + d.get("chkick", 0.0))
                ykick = bv * (d.get("vkick", 0.0) + d.get("cvkick", 0.0))
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if bool(d.get("sinkick", False)):
                    # TODO(audit): Map sinusoidally driven kicker modes (SINKICK/SINPEAK/SINTUNE/SINPHASE).
                    _warn(
                        f"KICKER '{d['name']}' uses sinusoidal kick options not translated; using static kick values."
                    )
                append_thin_with_length(
                    elements.Kicker(
                        name=d["name"],
                        xkick=xkick,
                        ykick=ykick,
                        rotation=tilt_degree,
                    ),
                    ds,
                    d["name"],
                )
            elif d["type"] == "hkicker":
                ds = d.get("l", 0.0)
                # MAD-X documents HKICKER with KICK, but its source also supports
                # HKICK as an alias and falls back to it when KICK is zero:
                # src/mad_dict.c (hkicker defines both kick and hkick) and
                # src/mad_elem.c:836-842 (el_par_value aliases kick <-> hkick).
                # Keep the alias support here so community MAD-X inputs such as
                # `HKICKER, HKICK=...` continue to translate correctly.
                # bv applies to the kicker kick (twiss/tracking).
                xkick0 = d.get("kick", 0.0)
                hkick0 = d.get("hkick", 0.0)
                if xkick0 == 0.0:
                    xkick0 = hkick0
                xkick = bv * (xkick0 + d.get("chkick", 0.0))
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if bool(d.get("sinkick", False)):
                    # TODO(audit): Map sinusoidally driven kicker modes (SINKICK/SINPEAK/SINTUNE/SINPHASE).
                    _warn(
                        f"HKICKER '{d['name']}' uses sinusoidal kick options not translated; using static kick values."
                    )
                append_thin_with_length(
                    elements.Kicker(
                        name=d["name"],
                        xkick=xkick,
                        ykick=0.0,
                        rotation=tilt_degree,
                    ),
                    ds,
                    d["name"],
                )
            elif d["type"] == "vkicker":
                ds = d.get("l", 0.0)
                # MAD-X likewise aliases VKICKER KICK <-> VKICK in el_par_value
                # (src/mad_elem.c:843-848). Preserve that behavior here so valid
                # MAD-X inputs using `VKICKER, VKICK=...` do not regress.
                ykick0 = d.get("kick", 0.0)
                vkick0 = d.get("vkick", 0.0)
                if ykick0 == 0.0:
                    ykick0 = vkick0
                ykick = bv * (ykick0 + d.get("cvkick", 0.0))
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if bool(d.get("sinkick", False)):
                    # TODO(audit): Map sinusoidally driven kicker modes (SINKICK/SINPEAK/SINTUNE/SINPHASE).
                    _warn(
                        f"VKICKER '{d['name']}' uses sinusoidal kick options not translated; using static kick values."
                    )
                append_thin_with_length(
                    elements.Kicker(
                        name=d["name"],
                        xkick=0.0,
                        ykick=ykick,
                        rotation=tilt_degree,
                    ),
                    ds,
                    d["name"],
                )
            elif d["type"] in ("monitor", "hmonitor", "vmonitor"):
                # MAD-X handles MONITOR/HMONITOR/VMONITOR with L>0 as a drift
                # (trrun.f90:820-821) and stores the reference orbit AFTER that
                # drift is applied (twiss.f90:1020 calls store_node_vector at
                # element exit). Emitting Drift(L) followed by a zero-length
                # ImpactX BeamMonitor therefore reproduces MAD-X's exit-face
                # recording convention exactly -- this is not a center-of-element
                # reading, despite the MAD-X `at=` convention pointing at the
                # element center.
                # TODO(audit): MAD-X MREX/MREY (monitor reading alignment errors
                # set via EALIGN) and HMONITOR/VMONITOR plane masking are not
                # translated. LRAD/TILT on MAD-X monitors are documented as
                # unused internally (changelog.tex:68), so we intentionally
                # ignore them here.
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
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    if any(abs(v) > 0.0 for v in (k2, k2s, tilt_degree)):
                        _warn(
                            f"SEXTUPOLE '{d['name']}' has L=0 and no integrated-strength form in MAD-X; using Marker."
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
                    impactx_beamline.append(elements.Marker(name=d["name"]))
                    if any(abs(v) > 0.0 for v in (k3, k3s, tilt_degree)):
                        _warn(
                            f"OCTUPOLE '{d['name']}' has L=0 and no integrated-strength form in MAD-X; using Marker."
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

                kick_components = []
                for order in range(max(max_order, -1) + 1):
                    kn = knl[order] if order < len(knl) else 0.0
                    ks = ksl[order] if order < len(ksl) else 0.0
                    if kn != 0.0 or ks != 0.0:
                        kick_components.append((order, kn, ks))

                if not kick_components:
                    # MULTIPOLE with no nonzero KNL/KSL: keep the element in the
                    # lattice as a Drift (if L>0) or Marker rather than silently
                    # dropping it.
                    if ds > 0.0:
                        impactx_beamline.append(
                            elements.Drift(name=d["name"], ds=ds, nslice=nslice)
                        )
                    else:
                        impactx_beamline.append(elements.Marker(name=d["name"]))
                else:
                    if ds > 0.0:
                        impactx_beamline.append(
                            elements.Drift(
                                name=d["name"] + "_drift_in",
                                ds=0.5 * ds,
                                nslice=nslice,
                            )
                        )
                    for order, kn, ks in kick_components:
                        impactx_beamline.append(
                            elements.Multipole(
                                name=f"{d['name']}_m{order + 1}",
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
                # MAD-X NLLENS lists `tilt` in the element dictionary
                # (mad_dict.c:3294) but neither tmnll (twiss.f90:8604) nor
                # ttnllens (trrun.f90:4014) reads it. Consume it here so we
                # don't emit a misleading "unread attribute" warning, and
                # explicitly note the MAD-X behavior if the user set one.
                if abs(d.get("tilt", 0.0)) > 0.0:
                    _warn(
                        f"NLLENS '{d['name']}' has nonzero TILT; MAD-X also ignores TILT for NLLENS (tmnll/ttnllens); not translated."
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
                # ImpactX ShortRF: V [1], freq [Hz], phase [deg]
                # ImpactX RFCavity: escale [1/m], freq [Hz], phase [deg]
                # MAD-X applies bvk to RFCAVITY VOLT (trrun.f90:1610, 1613:
                # rfv = bvk * node_value('volt ')).
                volt = bv * d.get("volt", 0.0)  # MV
                lag = d.get("lag", 0.0)  # fraction of 2pi
                lagtap = d.get("lagtap", 0.0)  # additional phase taper in turns
                # MAD-X dict default for HARMON is 0 (mad_dict.c:2034). The
                # effective frequency is then HARMON*FREQ0; both must be set
                # by the user for a usable cavity.
                harmon = d.get("harmon", 0.0)
                # MAD-X uses sin(2*pi*(lag - f*t)); ImpactX uses a cosine convention
                # with phase=0 on crest and phase=-90 at zero crossing.
                # TODO(audit): this matches MAD-X on crest, but the time-dependence
                # sign assumes ImpactX's t-coordinate has the same sign convention as
                # MAD-X's T (T = -c*dt). Verify with an off-crest bunching/debunching
                # round-trip against MAD-X PTC track before trusting longitudinal
                # focusing behavior. See ShortRF.H:160 and trrun.f90:1662.
                phase = (lag + lagtap) * 360.0 - 90.0
                tilt_degree = d.get("tilt", 0.0) * rad_to_deg
                if "freq" in d and d.get("freq", 0.0) > 0.0:
                    # MAD-X RFCAVITY frequency is specified in MHz.
                    freq = d.get("freq", 0.0) * 1.0e6
                else:
                    freq = harmon * freq0
                    if freq == 0.0:
                        # Corner case: FREQ unset and HARMON*FREQ0 yields 0
                        # (HARMON=0 by MAD-X default, or BEAM FREQ0=0).
                        # Keep translation deterministic (freq=0) and warn.
                        _warn(
                            f"RFCAVITY '{d['name']}' has no usable frequency: "
                            f"explicit FREQ unset and HARMON*FREQ0={harmon}*{freq0}=0; "
                            f"resulting RF frequency is 0."
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
                    if attr == "no_cavity_totalpath":
                        if attr in d and bool(d.get(attr, False)):
                            # MAD-X models this as a logical/PTC flag, not a numeric
                            # strength parameter.
                            _warn(
                                f"RFCAVITY '{d['name']}' attribute '{attr}' is not translated; using simplified model."
                            )
                    elif attr in d and abs(d.get(attr, 0.0)) > 0:
                        # TODO(audit): Extend translation for advanced MAD-X/PTC cavity attributes.
                        _warn(
                            f"RFCAVITY '{d['name']}' attribute '{attr}' is not translated; using simplified model."
                        )

                if ds > 0:
                    # Thick cavity: use RFCavity with a flat (pillbox) field profile.
                    # ImpactX evaluates the profile as 0.5*cos_coefficients[0] + ...
                    # (standard Fourier half-a0 convention, see RFCavity.H:385), so a
                    # uniform profile of 1 requires cos_coefficients=[2.0]. With that
                    # choice, the on-crest integrated energy gain is escale * ds, and
                    # escale = VOLT_norm / ds reproduces MAD-X's VOLT exactly.
                    escale = normalized_rf_escale(volt, ds)
                    impactx_beamline.append(
                        elements.RFCavity(
                            name=d["name"],
                            ds=ds,
                            escale=escale,
                            freq=freq,
                            phase=phase,
                            cos_coefficients=[2.0],
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
                            V=normalized_rf_voltage(volt),
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
                sorted(supported_madx_elements),
            )
    return impactx_beamline


def read_lattice(madx_file, nslice=1, *, line=None, sequence=None):
    """
    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param madx_file: file name to MAD-X file with beamline elements
    :param nslice: number of ds slices per element
    :param line: explicit MAD-X LINE name to expand when no USE command is present
    :param sequence: explicit MAD-X SEQUENCE/PERIOD name to expand
    :return: list of ImpactX.KnownElements
    """
    madx = MADXParser()
    madx.parse(madx_file)
    beamline = madx.getBeamline(line_name=line, sequence_name=sequence)
    freq0 = madx.getFreq0()
    if getattr(madx.context.beam, "mass", 0.0) > 0.0:
        ref_mass_MeV = madx.context.beam.mass * 1.0e3
    else:
        ref_mass_MeV = beam(
            madx.getParticle(),
            charge=getattr(madx.context.beam, "charge", None),
            mass=None,
            energy=float(madx.getEtot()),
        )["mass"]
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

    bv = float(getattr(madx.context.beam, "bv", 1.0))

    return lattice(
        beamline,
        nslice,
        freq0=freq0,
        options=options,
        ref_mass_MeV=ref_mass_MeV,
        bv=bv,
    )


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

    GeV2MeV = 1.0e3
    kg2MeV = sc.c**2 / sc.electron_volt * 1.0e-6
    electron_mass_MeV = sc.m_e * kg2MeV
    proton_mass_MeV = sc.m_p * kg2MeV
    muon_mass_MeV = (
        electron_mass_MeV / sc.physical_constants["electron-muon mass ratio"][0]
    )
    # Match MAD-X source (`mad_dict.c`), which defines ION with `mass = nmass`.
    ion_mass_MeV = 0.93956542052 * GeV2MeV
    if energy is None:
        energy_MeV = 1.0e3  # MAD-X default is 1 GeV total particle energy
    else:
        energy_MeV = energy * GeV2MeV
    if particle is None:
        particle = ""
    particle = particle.casefold()
    generic_mass_MeV = None if mass is None else mass * GeV2MeV

    impactx_beam = {
        "positron": {"mass": electron_mass_MeV, "charge": 1.0},
        "electron": {"mass": electron_mass_MeV, "charge": -1.0},
        "proton": {"mass": proton_mass_MeV, "charge": 1.0},
        "antiproton": {"mass": proton_mass_MeV, "charge": -1.0},
        "posmuon": {
            "mass": muon_mass_MeV,
            "charge": 1.0,
        },  # positively charged muon (anti-muon)
        "negmuon": {
            "mass": muon_mass_MeV,
            "charge": -1.0,
        },  # negatively charged muon
        "ion": {"mass": ion_mass_MeV, "charge": 1.0},
        "generic": {"mass": generic_mass_MeV, "charge": charge},
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
    if reference_particle["mass"] is None or reference_particle["charge"] is None:
        raise ValueError(
            "Unknown MAD-X particle species requires explicit MASS and CHARGE in the BEAM command."
        )
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
        charge=(
            None
            if getattr(madx.context.beam, "charge", 0.0) == 0.0
            else float(madx.context.beam.charge)
        ),
        mass=(
            None
            if getattr(madx.context.beam, "mass", 0.0) == 0.0
            else float(madx.context.beam.mass)
        ),
        energy=float(madx.getEtot()),  # MADX default energy is in GeV
    )

    ref.set_charge_qe(ref_particle_dict["charge"])
    ref.set_mass_MeV(ref_particle_dict["mass"])
    ref.set_kin_energy_MeV(ref_particle_dict["energy"] - ref_particle_dict["mass"])

    return ref
