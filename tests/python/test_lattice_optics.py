#!/usr/bin/env python3
#
# Copyright 2022-2024 The ImpactX Community
#
# Authors: Axel Huebl, Chad Mitchell
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from impactx import Config, Map6x6, RefPart, elements, push


def make_ref(*, x=0.0, y=0.0, px=0.0, py=0.0, species="electron", kin_energy_MeV=1.0e3):
    """Create a reference particle with configurable transverse coordinates."""
    ref = RefPart()
    ref.set_species(species).set_kin_energy_MeV(kin_energy_MeV)
    ref.x = x
    ref.y = y
    ref.px = px
    ref.py = py
    return ref


def map_equal(lhs, rhs):
    """Compare Map6x6 or NumPy-array-like values with precision-dependent tolerances."""
    lhs_np = lhs.to_numpy() if hasattr(lhs, "to_numpy") else lhs
    rhs_np = rhs.to_numpy() if hasattr(rhs, "to_numpy") else rhs

    if Config.precision == "SINGLE":
        return np.allclose(lhs_np, rhs_np, atol=1.0e-7, rtol=5.0e-5)
    else:
        return np.allclose(lhs_np, rhs_np, atol=1.0e-14, rtol=1.0e-8)


def sextupole_feeddown_quadrupole(ref, *, kn, ks, dx=0.0, dy=0.0, rotation=0.0):
    """Return the thin quadrupole equivalent to a sextupole linearized at ``ref``."""
    theta = np.deg2rad(rotation)
    xc = ref.x - dx
    yc = ref.y - dy
    zeta_local = complex(
        xc * np.cos(theta) + yc * np.sin(theta),
        -xc * np.sin(theta) + yc * np.cos(theta),
    )
    feeddown = complex(kn, ks) * zeta_local
    return elements.Multipole(
        multipole=2,
        K_normal=feeddown.real,
        K_skew=feeddown.imag,
        rotation=rotation,
    )


def test_lattice_linear_map():
    """Calculate the linear transfer map of a lattice."""

    # Create reference particle
    ref = make_ref()

    # Create a valid test lattice (all elements define a linear transfer map)
    lattice = elements.KnownElementsList(
        [
            elements.Sbend(name="bend1", ds=1.0, rc=10.0, nslice=2),
            elements.Drift(name="drift1", ds=2.0, nslice=2),
            elements.Quad(name="quad1", ds=0.5, k=1.0, nslice=2),
            elements.Drift(name="drift2", ds=1.0, nslice=2),
        ]
    )

    # Expected result (matrix multiplication)
    R_expected = np.array(
        [
            [3.74670546e-01, 2.54005827e00, 0, 0, 0, -2.34864805e-01],
            [-4.76219076e-01, -5.59489407e-01, 0, 0, 0, 3.20646253e-02],
            [0, 0, 1.64872127e00, 6.59488508e00, 0, 0],
            [0, 0, 5.21095305e-01, 2.69091188e00, 0, 0],
            [9.98334297e-02, 4.99583537e-02, 0, 0, 1.00000000e00, -1.66466013e-03],
            [0, 0, 0, 0, 0, 1.00000000e00],
        ]
    )

    # Calculate Linear Transfer Map
    R = lattice.transfer_map(ref)
    assert map_equal(R, R_expected)

    # Check unexpected/unsupported options
    with pytest.raises(RuntimeError):
        lattice.transfer_map(ref, order="invalid")


def test_multipole_on_axis_sextupole_has_identity_linear_map():
    """An on-axis thin sextupole has no linear term."""
    ref = make_ref()
    lattice = elements.KnownElementsList(
        [elements.Multipole(multipole=3, K_normal=2.5, K_skew=-0.7)]
    )

    assert map_equal(lattice.transfer_map(ref), Map6x6.identity())


def test_multipole_off_axis_sextupole_matches_feeddown_quadrupole():
    """A sextupole linearized about x!=0 must feed down to a thin quadrupole."""
    ref_offx = make_ref(x=2.0e-3)
    ref_center = make_ref()
    k2 = 11.0

    map_sextupole = elements.KnownElementsList(
        [elements.Multipole(multipole=3, K_normal=k2, K_skew=0.0)]
    ).transfer_map(ref_offx)

    map_quadrupole = elements.KnownElementsList(
        [elements.Multipole(multipole=2, K_normal=k2 * ref_offx.x, K_skew=0.0)]
    ).transfer_map(ref_center)

    assert map_equal(
        map_sextupole,
        map_quadrupole,
    )


def test_multipole_misaligned_sextupole_has_feeddown_at_zero_reference_offset():
    """A sextupole misalignment must generate feed-down even if ref.x=ref.y=0."""
    ref = make_ref()
    kn = 9.0
    ks = -4.0
    dx = 1.0e-3
    dy = -1.5e-3
    lattice_sextupole = elements.KnownElementsList(
        [elements.Multipole(multipole=3, K_normal=kn, K_skew=ks, dx=dx, dy=dy)]
    )

    lattice_quadrupole = elements.KnownElementsList(
        [sextupole_feeddown_quadrupole(ref, kn=kn, ks=ks, dx=dx, dy=dy)]
    )

    R = lattice_sextupole.transfer_map(ref)
    assert not map_equal(R, Map6x6.identity())
    assert map_equal(R, lattice_quadrupole.transfer_map(make_ref()))


def test_multipole_rolled_normal_quadrupole_matches_skew_quadrupole():
    """A normal quadrupole rolled by 45 degrees becomes a pure skew quadrupole."""
    ref = make_ref()
    k1 = 3.0

    map_rolled_normal = elements.KnownElementsList(
        [elements.Multipole(multipole=2, K_normal=k1, K_skew=0.0, rotation=45.0)]
    ).transfer_map(ref)

    map_skew = elements.KnownElementsList(
        [elements.Multipole(multipole=2, K_normal=0.0, K_skew=-k1)]
    ).transfer_map(ref)

    assert map_equal(map_rolled_normal, map_skew)


def test_multipole_reference_particle_ignores_alignment_errors():
    """Alignment errors act on the beam/map, not on the RefPart push itself."""
    ref_ideal = make_ref(x=2.0e-3, y=-1.0e-3, px=3.0e-4, py=-2.0e-4)
    ref_aligned = ref_ideal.copy()

    push(
        ref_ideal,
        elements.Multipole(multipole=3, K_normal=9.0, K_skew=-4.0),
    )
    push(
        ref_aligned,
        elements.Multipole(
            multipole=3,
            K_normal=9.0,
            K_skew=-4.0,
            dx=1.0e-3,
            dy=-1.5e-3,
            rotation=23.0,
        ),
    )

    for attr in ("x", "y", "px", "py", "pz", "t", "pt", "s"):
        assert np.isclose(getattr(ref_aligned, attr), getattr(ref_ideal, attr))


def test_multipole_feeddown_composes_across_multiple_misaligned_rotated_maps():
    """Misalignment/rotation feed-down must compose correctly across a lattice."""
    ref = make_ref(x=1.6e-3, y=-0.7e-3)

    drift = dict(ds=0.4, nslice=1)
    sextupole_1 = dict(kn=7.0, ks=-1.5, dx=0.9e-3, dy=-0.5e-3, rotation=17.0)
    sextupole_2 = dict(kn=-5.0, ks=2.25, dx=-1.1e-3, dy=0.8e-3, rotation=-31.0)

    lattice_sextupoles = elements.KnownElementsList(
        [
            elements.Multipole(
                multipole=3,
                K_normal=sextupole_1["kn"],
                K_skew=sextupole_1["ks"],
                dx=sextupole_1["dx"],
                dy=sextupole_1["dy"],
                rotation=sextupole_1["rotation"],
            ),
            elements.Drift(**drift),
            elements.Multipole(
                multipole=3,
                K_normal=sextupole_2["kn"],
                K_skew=sextupole_2["ks"],
                dx=sextupole_2["dx"],
                dy=sextupole_2["dy"],
                rotation=sextupole_2["rotation"],
            ),
        ]
    )

    lattice_feeddown_quads = elements.KnownElementsList(
        [
            sextupole_feeddown_quadrupole(ref, **sextupole_1),
            elements.Drift(**drift),
            sextupole_feeddown_quadrupole(ref, **sextupole_2),
        ]
    )

    assert map_equal(
        lattice_sextupoles.transfer_map(ref),
        lattice_feeddown_quads.transfer_map(ref),
    )
