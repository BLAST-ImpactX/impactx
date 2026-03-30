#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest


def test_default_contructor():
    io = synergia.collective.Impedance_options("foo.dat", "XLYLZ")


def test_constructor_onearg():
    io = synergia.collective.Impedance_options("foo.dat", "XLYLZ", 100)


def test_IO_defaults():
    io = synergia.collective.Impedance_options("foo.dat", "XLYLZ")

    assert io.z_grid == 1000

    assert not io.full_machine

    assert io.nstored_turns == 15

    assert io.num_buckets == 1

    assert io.orbit_length == 1.0

    assert io.bunch_spacing == 1.0

    assert io.mwf_xlead == 1.0

    assert io.mwf_xtrail == 1.0

    assert io.mwf_ylead == 1.0

    assert io.mwf_ytrail == 1.0

    assert io.mwf_zwake == 1.0


def test_IO_set_attributes():
    io = synergia.collective.Impedance_options("foo.dat", "XLYLZ", 100)
    assert io.z_grid == 100

    io.num_buckets = 4
    assert io.num_buckets == 4

    io.mwf_xlead = 0.5
    assert io.mwf_xlead == 0.5

    io.mwf_xtrail = 0.0
    assert io.mwf_xtrail == 0.0

    io.mwf_ylead = 0.25
    assert io.mwf_ylead == 0.25

    io.mwf_ytrail = -0.75
    assert io.mwf_ytrail == -0.75

    io.mwf_zwake = 0.375
    assert io.mwf_zwake == 0.375
