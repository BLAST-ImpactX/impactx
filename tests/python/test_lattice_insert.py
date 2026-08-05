#!/usr/bin/env python3
#
# Copyright 2022-2026 The ImpactX Community
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-


from impactx import elements


def test_element_insert():
    """
    This tests the lattice element insert every ds elements.
    """
    monitor = elements.BeamMonitor("monitor", backend="h5")

    fodo = [
        elements.Drift(name="drift1", ds=0.25),
        elements.Quad(name="quad1", ds=1.0, k=1.0),
        elements.Drift(name="drift2", ds=0.5),
        elements.Quad(name="quad2", ds=1.0, k=-1.0),
        elements.Drift(name="drift3", ds=0.25),
    ]

    fodo = elements.transformation.insert_element_every_ds(
        # todo: avoid that we need to cast
        elements.KnownElementsList(fodo),
        ds=0.1,
        element=monitor,
    )

    inserted_monitors = list(
        filter(lambda el: el.to_dict()["type"] == "BeamMonitor", fodo)
    )

    assert len(fodo) == 63
    assert len(inserted_monitors) == 29

    # clean shutdown
    monitor.finalize()


def test_element_insert_unnamed():
    """Splitting an unnamed element gives its leftover segments generated names."""
    lattice = elements.transformation.insert_element_every_ds(
        elements.KnownElementsList([elements.Drift(ds=0.25)]),
        ds=0.1,
        element=elements.Empty(),
    )

    drifts = [el for el in lattice if el.to_dict()["type"] == "Drift"]

    assert len(lattice) == 5
    assert [drift.name for drift in drifts] == [
        None,
        "_leftover",
        "_leftover_leftover",
    ]
