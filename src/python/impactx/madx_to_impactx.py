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


def lattice(parsed_beamline, nslice=1, freq0=0.0):
    """
    Function that converts a list of elements in the MADXParser format into ImpactX format
    :param parsed_beamline: list of dictionaries
    :param nslice: number of ds slices per element
    :param freq0: revolution frequency in MHz (from BEAM command)
    :return: list of translated dictionaries
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
        "MULTIPOLE": "Multipole",
        "SEXTUPOLE": "ExactMultipole",
        "OCTUPOLE": "ExactMultipole",
        "RFCAVITY": "RFCavity",
        "NLLENS": "NonlinearLens",
        # TODO Figure out how to identify these
        "ShortRF": "ShortRF",
        "ConstF": "ConstF",
    }

    impactx_beamline = []

    for d in parsed_beamline:
        # print(d)
        if d["type"] in [k.casefold() for k in list(madx_to_impactx_dict.keys())]:
            if d["type"] == "drift":
                impactx_beamline.append(
                    elements.Drift(name=d["name"], ds=d["l"], nslice=nslice)
                )
            elif d["type"] == "quadrupole":
                impactx_beamline.append(
                    elements.Quad(name=d["name"], ds=d["l"], k=d["k1"], nslice=nslice)
                )
            elif d["type"] == "sbend":
                impactx_beamline.append(
                    elements.Sbend(
                        name=d["name"], ds=d["l"], rc=d["l"] / d["angle"], nslice=nslice
                    )
                )
            elif d["type"] == "rbend":
                # RBEND: rectangular bending magnet
                # Converts to DipEdge + SBend + DipEdge
                # Reference: https://mad.web.cern.ch/mad/webguide/manual.html#Ch11.S3
                angle = d["angle"]
                l_chord = d["l"]

                # Convert chord length to arc length (RBARC=true is MAD-X default)
                if abs(angle) > 0:
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
                    )
                )
                # SBend body
                impactx_beamline.append(
                    elements.Sbend(
                        name=d["name"],
                        ds=l_arc,
                        rc=rc,
                        nslice=nslice,
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
                    )
                )
            elif d["type"] == "solenoid":
                impactx_beamline.append(
                    elements.Sol(name=d["name"], ds=d["l"], ks=d["ks"], nslice=nslice)
                )
            elif d["type"] == "dipedge":
                h = d.get("h", 0.0)
                impactx_beamline.append(
                    elements.DipEdge(
                        name=d["name"],
                        psi=d.get("e1", 0.0),
                        rc=1.0 / h if h != 0.0 else 0.0,
                        # MAD-X is using half the gap height
                        g=2.0 * d.get("hgap", 0.0),
                        K2=d.get("fint", 0.0),
                    )
                )
            elif d["type"] in ("kicker", "tkicker"):
                impactx_beamline.append(
                    elements.Kicker(
                        name=d["name"],
                        xkick=d.get("hkick", 0.0),
                        ykick=d.get("vkick", 0.0),
                    )
                )
            elif d["type"] == "hkicker":
                impactx_beamline.append(
                    elements.Kicker(
                        name=d["name"],
                        xkick=d.get("hkick", 0.0),
                        ykick=0.0,
                    )
                )
            elif d["type"] == "vkicker":
                impactx_beamline.append(
                    elements.Kicker(
                        name=d["name"],
                        xkick=0.0,
                        ykick=d.get("vkick", 0.0),
                    )
                )
            elif d["type"] == "monitor":
                if d["l"] > 0:
                    impactx_beamline.append(
                        elements.Drift(
                            name=d["name"] + "_drift", ds=d["l"], nslice=nslice
                        )
                    )
                impactx_beamline.append(
                    elements.BeamMonitor(name="monitor", backend="h5")
                )  # TODO: use name=d["name"] ?
            elif d["type"] == "sextupole":
                k2 = d.get("k2", 0.0)
                impactx_beamline.append(
                    elements.ExactMultipole(
                        name=d["name"],
                        ds=d["l"],
                        k_normal=[0.0, 0.0, k2],
                        k_skew=[0.0, 0.0, 0.0],
                        nslice=nslice,
                    )
                )
            elif d["type"] == "octupole":
                k3 = d.get("k3", 0.0)
                impactx_beamline.append(
                    elements.ExactMultipole(
                        name=d["name"],
                        ds=d["l"],
                        k_normal=[0.0, 0.0, 0.0, k3],
                        k_skew=[0.0, 0.0, 0.0, 0.0],
                        nslice=nslice,
                    )
                )
            elif d["type"] == "rfcavity":
                # MAD-X: volt [MV], lag [2pi], harmon [1]
                # ImpactX ShortRF: V [MV], freq [Hz], phase [deg]
                # ImpactX RFCavity: escale [MV/m], freq [Hz], phase [deg]
                volt = d.get("volt", 0.0)  # MV
                lag = d.get("lag", 0.0)  # fraction of 2pi
                harmon = d.get("harmon", 1.0)
                phase = lag * 360.0  # convert to degrees
                freq = harmon * freq0 * 1.0e6  # MHz -> Hz
                ds = d.get("l", 0.0)

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
                        )
                    )
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

    return lattice(beamline, nslice, freq0=freq0)


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
        warnings.warn(
            f"Particle species name '{particle}' not in '{impactx_beam.keys()}'",
            UserWarning,
        )
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
    :param RefPart ref: ImpactX reference particle (passed by reference)
    :param madx_file: file name to MAD-X file with beamline elements
    :return: list of ImpactX.KnownElements
    """
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
